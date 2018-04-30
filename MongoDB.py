import sys
import itertools
try:
    from itertools import izip_longest as zip_longest
except ImportError:
    from itertools import zip_longest
import os
from datetime import datetime
from io import BytesIO as IOFunc
import zstd
from collections import defaultdict

from Bio import SeqIO

import asyncio
import motor.motor_asyncio
import pymongo

from SwissProtUtils import parse_raw_swiss, filter_proks
from tqdm import tqdm


def _get_date(dateline):
    months = {
        'JAN': 1, 'FEB': 2, 'MAR': 3,
        'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9,
        'OCT': 10, 'NOV': 11, 'DEC': 12,
    }
    day, month, year = dateline.split()[1].strip(',').split('-')
    return datetime(int(year), months[month], int(day))

def _create_protein(raw_record, compressor):
    lines = raw_record.decode().split('\n')
    desc_lines = []
    refs = defaultdict(list)
    for l in lines:
        s = l[:2]
        if s == 'CC' or s == '  ':
            continue
        elif s == 'DT':
            dateline = l
        elif s == 'DE':
            desc_lines.append(l.split(maxsplit=1)[1])
        elif s == 'OS':
            genome = l.split(maxsplit=1)[1].strip('. ')
        elif s == 'OX':
            taxid = l.split('=')[1].strip(';')
        elif s == 'DR':
            ref = l.split(maxsplit=1)[1]
            dec = ref.strip('.').split(';')
            refs[dec[0]].append(dec[1].strip())

    return dict(
        _id=lines[1].split()[1].strip(';'),
        genome=genome,
        taxid=taxid,
        description=' '.join(desc_lines),
        updated=_get_date(dateline),
        raw_record=compressor.compress(raw_record),
        Uni_name=[lines[0].split()[1]],
        **refs,
    )


class MongoDatabase(object):

    ids = ['_id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name']

    def __init__(self, database, host):
        self.loop = asyncio.get_event_loop()
        self.host = host
        self.client = motor.motor_asyncio.AsyncIOMotorClient(*host)
        self.database = database
        self.col = self.client[database].proteins
        if self.loop.run_until_complete(self.col.count({'_id': 'dict_data'})):
            dict_data = self.loop.run_until_complete(self.client[database].proteins.find_one('dict_data'))['dict_data']
            self.compressor = zstd.ZstdCompressor(dict_data=zstd.ZstdCompressionDict(dict_data))
            self.decomp = zstd.ZstdDecompressor(dict_data=zstd.ZstdCompressionDict(dict_data))
        else:
            self.initialize([], database=self.database)

    def get_item(self, item):
        t = self.loop.run_until_complete(self.col.find_one({'$or': [{i: item} for i in self.ids]}))
        if t is None:
            return None
        r = SeqIO.read(IOFunc(self.decomp.decompress(t['raw_record'])), 'swiss')
        return r


    def get_iter(self):
        i = self.loop.run_until_complete(self._get_iter())
        while not i.empty():
            yield self.loop.run_until_complete(i.get())

    async def _get_iter(self):
        q = asyncio.Queue()
        async for entry in self.col.find({'Uni_name': {'$exists': True}}):
            await q.put(SeqIO.read(IOFunc(self.decomp.decompress(entry['raw_record'])), 'swiss'))
        return q

    def get_iterkeys(self):
        i = self.loop.run_until_complete(self._get_iterkeys())
        while not i.empty():
            yield self.loop.run_until_complete(i.get())

    async def _get_iterkeys(self):
        q = asyncio.Queue()
        async for i in self.col.find({'Uni_name': {'$exists': True}}, {'_id': 1}):
            await q.put(i['_id'])
        return q

    def get_keys(self):
        return self.loop.run_until_complete(self.col.distinct('_id'))


    def length(self):
        return self.loop.run_until_complete(self.col.count({'Uni_name': {'$exists': True}}))


    def get_by(self, attr, value):
        return self.loop.run_until_complete(self._get_by(attr, value))

    async def _get_by(self, attr, value):
        ret = []
        res = self.col.find({attr: value}, {'raw_record': 1})
        async for i in res:
            ret.append(SeqIO.read(IOFunc(self.decomp.decompress(i['raw_record'])), 'swiss'))
        return ret


    def initialize(self, seq_handles, database='uniprot', filter_fn=None, loud=False, n_seqs=None):
        if loud:
            print("--initializating database\n", file=sys.stderr)
        self.loop.run_until_complete(self.client[database].proteins.drop())
        with open(os.path.dirname(__file__)+"/cc16.zstd_dict", 'rb') as i:
            d = i.read()
            self.loop.run_until_complete(self.client[database].proteins.insert_one({'_id': 'dict_data', 'dict_data': d}))

        self.loop.run_until_complete(self.add_from_handles(seq_handles, filter_fn))

        self.loop.run_until_complete(self.client[database].proteins.create_index([('genome', 1)]))
        indices = ['RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'PDB', 'PDBsum', 'EMBL', 'KEGG',
                   'GO', 'InterPro', 'Pfam', 'TIGRFAMs', 'PROSITE']
        for field in indices:
            self.loop.run_until_complete(self.client[database].proteins.create_index(keys=[(field, pymongo.ASCENDING)],
                                                   partialFilterExpression={field: {'$exists': True}}))
        if loud:
            print("--initialized database\n", file=sys.stderr)

    def add_protein(self, raw_record, test=None, test_attr=None):
        protein = _create_protein(raw_record, self.compressor)
        self.loop.run_until_complete(self._add_protein(protein, test=None, test_attr=None))

    async def _add_protein(self, protein, test=None, test_attr=None):
        if test:
            good = False
            if test == protein['_id']:
                good = True
            for ref in [test_attr] if test_attr else self.ids:
                if test in getattr(protein, ref, []):
                    good = True
            if not good:
                return False

        await self.col.replace_one({'_id': protein['_id']}, protein, upsert=True)

    def update(self, handles, filter_fn=None, loud=False, total=None):
        self.loop.run_until_complete(self.add_from_handles(handles, filter_fn=filter_fn, total=total))


    async def add_from_handles(self, handles, filter_fn=None, total=None):
        raw_protein_records = itertools.chain(*[parse_raw_swiss(handle, filter_fn, check_date=True) for handle in handles])
        additions = []
        with tqdm(total=total, smoothing=0.1) as pbar:
            for record, date in raw_protein_records:
                #if date > getattr(await self.col.find_one(record.split(maxsplit=2)[1]), 'updated', datetime.min):
                protein = _create_protein(record, self.compressor)
                additions.append(self.loop.create_task(self._add_protein(protein)))
                pbar.update()
        return asyncio.gather(*additions)
