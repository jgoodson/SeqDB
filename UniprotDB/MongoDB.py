import itertools
from typing import Union, Iterable, Callable, Generator, List

import pymongo
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from UniprotDB.BaseDatabase import BaseDatabase
from UniprotDB.SwissProtUtils import parse_raw_swiss


class MongoDatabase(BaseDatabase):

    def __init__(self, database: str, host: Union[str, tuple] = ('localhost',), **kwargs):
        super().__init__(database, host, **kwargs)
        self.client = pymongo.MongoClient(*host)
        self.col = self.client[database].proteins

    def get_item(self, item: str) -> Union[SeqRecord, None]:
        t = self.col.find_one({'$or': [{i: item} for i in self.ids]}, {'raw_record': True})
        if t is None:
            return None
        r = self._extract_seqrecord(t['raw_record'])
        return r

    def get_iter(self) -> Generator[SeqRecord, None, None]:
        for entry in self.col.find({}, {'raw_record': True}):
            yield self._extract_seqrecord(entry['raw_record'])

    def get_iterkeys(self) -> Generator[str, None, None]:
        for i in self.col.find({}, {'_id': True}):
            yield i['_id']

    def get_keys(self) -> List[str]:
        return self.col.distinct('_id')

    def length(self) -> int:
        return self.col.count_documents({})

    def get_by(self, attr: str, value: str) -> List[SeqRecord]:
        ret = []
        res = self.col.find({attr: value}, {'raw_record': True})
        for i in res:
            ret.append(self._extract_seqrecord(i['raw_record']))
        return ret

    def _reset(self) -> None:
        self.client[self.database].proteins.drop()

    def _create_indices(self, background: bool = False) -> None:
        for field in self.indices:
            self.client[self.database].proteins.create_index([(field, pymongo.ASCENDING)],
                                                             background=background, sparse=True)

    def update(self, handles: Iterable, filter_fn: Callable = None,
               loud: bool = False, total: int = None, workers: int = 1, fake: bool = False) -> None:
        raw_protein_records = itertools.chain(*[parse_raw_swiss(handle, filter_fn) for handle in handles])
        for record in tqdm(raw_protein_records, disable=(not loud), total=total, smoothing=0.1):
            if not fake:
                self.add_protein(self.create_protein_func(record))

    def add_protein(self, protein: dict) -> bool:
        self.col.replace_one({'_id': protein['_id']}, protein, upsert=True)
        return True
