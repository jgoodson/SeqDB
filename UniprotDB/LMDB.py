import itertools

import lmdb
import bson
from tqdm import tqdm

from UniprotDB.BaseDatabase import BaseDatabase
from UniprotDB.SwissProtUtils import parse_raw_swiss
from UniprotDB._utils import _create_record_swiss

class LMDBDatabase(BaseDatabase):

    def __init__(self, database, host='seqdb.lmdb', index=True, map_size=int(1024**4), **kwargs):
        super().__init__(database, host, **kwargs)
        self.client = lmdb.open(host, max_dbs=1+len(self.indices), map_size=map_size)
        self.has_index = index
        self._setup_dbs()


    def _setup_dbs(self):
        self.db = self.client.open_db(self.database.encode())
        if self.has_index:
            self.index_dbs = {}
            for index in self.indices:
                self.index_dbs[index] = self.client.open_db(index.encode(), dupsort=True)

    def _reset(self):
        with self.client.begin(write=True) as txn:
            txn.drop(self.db)
            if self.has_index:
                for index in self.indices:
                    txn.drop(self.index_dbs[index])
        self._setup_dbs()

    def get_item(self, item):
        with self.client.begin() as txn:
            t = txn.get(item.encode(), db=self.db)
        if t is None:
            return None
        return self._extract_seqrecord(t)

    def get_iter(self):
        with self.client.begin() as txn:
            cursor = txn.cursor(self.db)
            for entry in cursor.iternext(keys=False):
                yield self._extract_seqrecord(entry)

    def get_iterkeys(self):
        with self.client.begin() as txn:
            cursor = txn.cursor(self.db)
            for entry in cursor.iternext(values=False):
                yield entry.decode()

    def get_keys(self):
        return list(self.get_iterkeys())

    def length(self):
        with self.client.begin() as txn:
            return txn.stat(self.db)['entries']

    def get_by(self, attr, value):
        ret = []
        with self.client.begin() as txn:
            if attr == '_id':
                ret.append(self.get_item(value))
                return ret
            if self.has_index:
                cur = txn.cursor(db=self.index_dbs[attr])
                if cur.set_key(value.encode()):
                    for i in cur.iternext_dup():
                        ret.append(self.get_item(i.decode()))
        return ret

    def _create_indices(self, background=False):
        pass

    def update(self, handles, filter_fn=None, loud=False, total=None, workers=1):
        self._add_from_handles(handles, filter_fn=filter_fn, total=total, loud=loud)

    def add_protein(self, protein):
        pid = protein['_id'].encode()
        with self.client.begin(write=True) as txn:
            txn.put(pid, protein['raw_record'], db=self.db)
            if self.has_index:
                for index in self.indices:
                    if index in protein:
                        if isinstance(protein[index], list):
                            for idx in protein[index]:
                                txn.put(idx.encode(), pid, db=self.index_dbs[index])
                        elif isinstance(protein[index], str):
                            txn.put(protein[index].encode(), pid, db=self.index_dbs[index])
                        else:
                            txn.put(str(protein[index]).encode(), pid, db=self.index_dbs[index])
        return True

    def _add_from_handles(self, handles, filter_fn=None, total=None, loud=False, fake=False):
        raw_protein_records = itertools.chain(*[parse_raw_swiss(handle, filter_fn) for handle in handles])
        for record in tqdm(raw_protein_records, disable=(not loud), total=total, smoothing=0.1):
            if not fake:
                self.add_protein(self.create_protein_func(record))
