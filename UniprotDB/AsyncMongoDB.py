import asyncio
import itertools

import motor.motor_asyncio
import pymongo
from tqdm import tqdm

from UniprotDB.BaseDatabase import BaseDatabase
from UniprotDB.SwissProtUtils import parse_raw_swiss


class MongoDatabase(BaseDatabase):

    def __init__(self, database, host, **kwargs):
        super().__init__(database, host, **kwargs)
        try:
            self.loop = asyncio.get_event_loop()
        except:
            self.loop = asyncio.new_event_loop()
            asyncio.set_event_loop(self.loop)
        self.client = motor.motor_asyncio.AsyncIOMotorClient(*host)
        self.col = self.client[database].proteins

    def get_item(self, item):
        t = self.loop.run_until_complete(self.col.find_one({'$or': [{i: item} for i in self.ids]}))
        if t is None:
            return None
        r = self._extract_seqrecord(t['raw_record'])
        return r

    def get_iter(self):
        q = asyncio.Queue()
        self.loop.create_task(self._get_iter(q))
        r = self.loop.run_until_complete(q.get())
        while r:
            yield r
            r = self.loop.run_until_complete(q.get())

    async def _get_iter(self, q):
        async for entry in self.col.find({'_id': {'$exists': True}}):
            await q.put(self._extract_seqrecord(entry['raw_record']))
        await q.put(None)

    def get_iterkeys(self):
        q = asyncio.Queue()
        self.loop.create_task(self._get_iterkeys(q))
        r = self.loop.run_until_complete(q.get())
        while r:
            yield r
            r = self.loop.run_until_complete(q.get())

    async def _get_iterkeys(self, q):
        async for i in self.col.find({'_id': {'$exists': True}}, {'_id': 1}):
            await q.put(i['_id'])
        await q.put(None)

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
            ret.append(self._extract_seqrecord(i['raw_record']))
        return ret

    def _reset(self):
        self.loop.run_until_complete(self.client[self.database].proteins.drop())

    def _create_indices(self, background=False):
        for field in self.indices:
            self.loop.run_until_complete(
                self.client[self.database].proteins.create_index([(field, pymongo.ASCENDING)],
                                                                 background=background, sparse=True))

    def update(self, handles, filter_fn=None, loud=False, total=None, workers=1):
        self.loop.run_until_complete(
            self._add_from_handles(handles, filter_fn=filter_fn, total=total, loud=loud)
        )

    def add_protein(self, protein):
        return self.loop.run_until_complete(self._add_protein(protein))

    async def _add_protein(self, protein):
        await self.col.replace_one({'_id': protein['_id']}, protein, upsert=True)
        return True

    async def _add_from_handles(self, handles, filter_fn=None, total=None, loud=False, fake=False):
        raw_protein_records = itertools.chain(*[parse_raw_swiss(handle, filter_fn) for handle in handles])
        for record in tqdm(raw_protein_records, disable=(not loud), total=total, smoothing=0.1):
            if not fake:
                await self._add_protein(self.create_protein_func(record))