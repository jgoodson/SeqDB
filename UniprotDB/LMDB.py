import itertools
import os.path
from typing import Iterable, Callable, Generator, List, BinaryIO, Union

import lmdb
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from UniprotDB.BaseDatabase import BaseDatabase
from UniprotDB.SwissProtUtils import parse_raw_swiss
from UniprotDB._utils import _create_protein_swiss_bytes


class RawLMDBDatabase(BaseDatabase):
    # TODO: This exhibits huge write slowdown. We split the main db into ten by last-digit which helped some.
    def __init__(self, database: str,
                 host: str = '~/.seqdb/',
                 index: bool = True,
                 map_size: int = int(1024 ** 4),
                 **kwargs):
        if host.startswith('~'):
            host = os.path.expanduser(host)
        try:
            os.mkdir(host)
        except FileExistsError:
            pass
        self.host = host
        self.map_size = map_size
        super().__init__(database, host, create_protein_func=_create_protein_swiss_bytes, **kwargs)
        self.has_index = index
        self._setup_dbs()

    def _setup_dbs(self) -> None:
        self.db = {}
        for i in range(10):
            self.db[str(i).encode()] = lmdb.open(os.path.join(self.host, str(i) + '.lmdb'), map_size=self.map_size,
                                                 writemap=True, map_async=True, readahead=False)
        if self.has_index:
            self.index_dbs = {}
            for index in self.indices:
                self.index_dbs[index] = lmdb.open(os.path.join(self.host, index + '.lmdb'), map_size=self.map_size,
                                                  writemap=True, map_async=True, readahead=False)

    def _reset(self) -> None:
        for i in range(10):
            with self.db[str(i).encode()].begin(write=True) as txn:
                db = self.db[str(i).encode()].open_db()
                txn.drop(db)
        if self.has_index:
            for index in self.indices:
                with self.index_dbs[index].begin(write=True) as txn:
                    db = self.index_dbs[index].open_db()
                    txn.drop(db)
        self._setup_dbs()

    def get_item(self, item: str) -> Union[SeqRecord, None]:
        with self.db[item[-1:].encode()].begin() as txn:
            t = txn.get(item.encode())
        if t is None:
            return None
        return self._extract_seqrecord(t)

    def get_iter(self) -> Generator[SeqRecord, None, None]:
        for i in range(10):
            with self.db[str(i).encode()].begin() as txn:
                cursor = txn.cursor()
                for entry in cursor.iternext(keys=False):
                    yield self._extract_seqrecord(entry)

    def get_iterkeys(self) -> Generator[str, None, None]:
        for i in range(10):
            with self.db[str(i).encode()].begin() as txn:
                cursor = txn.cursor()
                for entry in cursor.iternext(values=False):
                    yield entry.decode()

    def get_keys(self) -> List[str]:
        return list(self.get_iterkeys())

    def length(self) -> int:
        total = 0
        for i in range(10):
            with self.db[str(i).encode()].begin() as txn:
                total += txn.stat()['entries']
        return total

    def get_by(self, attr: str, value: str) -> List[SeqRecord]:
        ret = []
        if attr == '_id':
            ret.append(self.get_item(value))
            return ret
        if self.has_index:
            with self.index_dbs[attr].begin() as txn:
                db = self.index_dbs[attr].open_db(dupsort=True)
                cur = txn.cursor(db=db)
                if cur.set_key(value.encode()):
                    for i in cur.iternext_dup():
                        ret.append(self.get_item(i.decode()))
        return ret

    def _create_indices(self, background: bool = False) -> None:
        pass

    def update(self, handles: Iterable, filter_fn: Callable = None,
               loud: bool = False, total: int = None, workers: int = 1) -> None:
        self._add_from_handles(handles, filter_fn=filter_fn, total=total, loud=loud)

    def add_protein(self, protein: dict) -> bool:
        pid = protein['_id']
        with self.db[pid[-1:]].begin(write=True) as txn:
            txn.put(pid, protein['raw_record'])
        if self.has_index:
            for index in self.indices:
                if index in protein:
                    with self.index_dbs[index].begin(write=True) as txn:
                        db = self.index_dbs[index].open_db(dupsort=True)
                        if isinstance(protein[index], list):
                            for idx in protein[index]:
                                txn.put(idx, pid, db=db)
                        elif isinstance(protein[index], bytes):
                            txn.put(protein[index], pid, db=db)
                        else:
                            txn.put(str(protein[index]).encode(), pid, db=db)
        return True

    def _add_from_handles(self, handles: Iterable[BinaryIO], filter_fn: Callable = None,
                          total: int = None, loud: bool = False, fake: bool = False) -> None:
        raw_protein_records = itertools.chain(*[parse_raw_swiss(handle, filter_fn) for handle in handles])
        for record in tqdm(raw_protein_records, disable=(not loud), total=total, smoothing=0.1):
            if not fake:
                self.add_protein(self.create_protein_func(record))


class LMDBDatabase(BaseDatabase):

    def __init__(self, database: str,
                 host: str,
                 **kwargs):
        import requests
        self.r = requests
        super().__init__(database, host, create_protein_func=_create_protein_swiss_bytes, **kwargs)
        self.host = host
        self.database = database

    def get_item(self, item: str) -> Union[SeqRecord, None]:
        t = self.r.get(self.host + '/get_item', params={'item': item}).json()
        if t['item'] is '':
            return None
        return self._extract_seqrecord(t[item])

    def get_iter(self) -> Generator[SeqRecord, None, None]:
        return NotImplementedError

    def get_iterkeys(self) -> Generator[str, None, None]:
        return NotImplementedError

    def get_keys(self) -> List[str]:
        return NotImplementedError

    def length(self) -> int:
        t = self.r.get(self.host + '/get_length').json()
        return t['length']

    def get_by(self, attr: str, value: str) -> List[SeqRecord]:
        t = self.r.get(self.host + '/get_by', params={'attr': attr, 'value': value}).json()
        ret = t['values']
        return ret

    def _create_indices(self, background: bool = False) -> None:
        pass

    def update(self, handles: Iterable, filter_fn: Callable = None,
               loud: bool = False, total: int = None, workers: int = 1) -> None:
        self._add_from_handles(handles, filter_fn=filter_fn, total=total, loud=loud)

    def add_protein(self, protein: bytes) -> bool:
        r = self.r.put(self.host + '/put_protein', params={'raw_record': protein})
        return r.status_code in (200, 204)

    def _add_from_handles(self, handles: Iterable[BinaryIO], filter_fn: Callable = None,
                          total: int = None, loud: bool = False, fake: bool = False) -> None:
        raw_protein_records = itertools.chain(*[parse_raw_swiss(handle, filter_fn) for handle in handles])
        for record in tqdm(raw_protein_records, disable=(not loud), total=total, smoothing=0.1):
            if not fake:
                self.add_protein(record)
