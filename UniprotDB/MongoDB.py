import itertools
import sys

try:
    from itertools import izip_longest as zip_longest
except ImportError:
    from itertools import zip_longest

from io import StringIO as IOFunc
import zstd

from Bio import SeqIO

import concurrent.futures
import asyncio
import motor.motor_asyncio
import pymongo
from tqdm import tqdm

from UniprotDB.SwissProtUtils import parse_raw_swiss
from UniprotDB._utils import _extract_seqrecord, _create_protein, _get_date

class MongoDatabase(object):

    ids = ['_id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name']

    def __init__(self, database, host):
        self.loop = asyncio.get_event_loop()
        #self.loop.set_debug(True)
        self.host = host
        self.client = motor.motor_asyncio.AsyncIOMotorClient(*host)
        self.database = database
        self.col = self.client[database].proteins

    def get_item(self, item):
        t = self.loop.run_until_complete(self.col.find_one({'$or': [{i: item} for i in self.ids]}))
        if t is None:
            return None
        r = _extract_seqrecord(t['raw_record'])
        return r


    def get_iter(self):
        i = self.loop.run_until_complete(self._get_iter())
        while not i.empty():
            yield self.loop.run_until_complete(i.get())

    async def _get_iter(self):
        async for entry in self.col.find({'_id': {'$exists': True}}):
            yield _extract_seqrecord(entry['raw_record'])

    def get_iterkeys(self):
        i = self.loop.run_until_complete(self._get_iterkeys())
        while not i.empty():
            yield self.loop.run_until_complete(i.get())

    async def _get_iterkeys(self):
        q = asyncio.Queue()
        async for i in self.col.find({'_id': {'$exists': True}}, {'_id': 1}):
            await q.put(i['_id'])
        return q

    def get_keys(self):
        return self.loop.run_until_complete(self.col.distinct('_id'))


    def length(self):
        return self.loop.run_until_complete(self.col.count_documents({'_id': {'$exists': True}}))

    def get_by(self, attr, value):
        return self.loop.run_until_complete(self._get_by(attr, value))

    async def _get_by(self, attr, value):
        ret = []
        res = self.col.find({attr: value}, {'raw_record': 1})
        async for i in res:
            ret.append(_extract_seqrecord(i['raw_record']))
        return ret


    def initialize(self, seq_handles, filter_fn=None, loud=False, n_seqs=None, processes=1):
        if loud:
            print("--initializating database\n", file=sys.stderr)
        self.loop.run_until_complete(self.client[self.database].proteins.drop())

        self.loop.run_until_complete(self.add_from_handles(seq_handles, filter_fn=filter_fn, total=n_seqs, processes=processes))

        self.loop.run_until_complete(self.client[self.database].proteins.create_index([('genome', 1)]))
        indices = ['RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'PDB', 'EMBL', 'GO', 'Pfam', 'Proteomes']
        for field in indices:
            self.loop.run_until_complete(self.client[self.database].proteins.create_index(keys=[(field, pymongo.ASCENDING)],
                                                   partialFilterExpression={field: {'$exists': True}}))
        if loud:
            print("--initialized database\n", file=sys.stderr)

    def add_protein(self, raw_record, test=None, test_attr=None):
        return self.loop.run_until_complete(self._add_protein(raw_record, test=test, test_attr=test_attr))

    async def _add_protein(self, record, ppe=None, test=None, test_attr=None):
        if ppe:
            protein = await self.loop.run_in_executor(ppe, _create_protein, record)
        else:
            protein = _create_protein(record)
        if test:
            good = False
            if test == protein['_id']:
                good = True
            if not good:
                for ref in ([test_attr] if test_attr else self.ids):
                    if test in protein.get(ref, []):
                        good = True
            if not good:
                return False
        await self.col.replace_one({'_id': protein['_id']}, protein, upsert=True)
        return True

    def update(self, handles, filter_fn=None, loud=False, total=None, processes=1):
        self.loop.run_until_complete(self.add_from_handles(handles, filter_fn=filter_fn, total=total, loud=loud, processes=processes))


    async def add_from_handles(self, handles, filter_fn=None, total=None, loud=False, processes=4):
        raw_protein_records = itertools.chain(*[parse_raw_swiss(handle, filter_fn) for handle in handles])
        tasks = []
        n = 100
        ppe = concurrent.futures.ProcessPoolExecutor(max_workers=processes)
        with tqdm(total=total, smoothing=0.1, disable=(not loud)) as pbar:
            for record in raw_protein_records:
                if len(tasks)>n:
                    done, pending = await asyncio.wait(tasks, return_when=asyncio.FIRST_COMPLETED)
                    for d in done:
                        tasks.remove(d)
                    if len(pending)>n:
                        for i in range(len(pending)-n):
                            await pending[i]
                tasks.append(asyncio.ensure_future(self._add_protein(record, ppe)))
                pbar.update()
        await asyncio.wait(tasks)
