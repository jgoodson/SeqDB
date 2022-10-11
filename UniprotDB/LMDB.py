import hashlib
import itertools
import json
import os.path
import shutil
from typing import Iterable, Callable, Generator, List, BinaryIO, Union, Dict

import lmdb
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from UniprotDB.BaseDatabase import BaseDatabase
from UniprotDB.SwissProtUtils import parse_raw_swiss


class RawLMDBDatabase(BaseDatabase):
    def __init__(self, database: str,
                 host: str = '~/.seqdb/',
                 index: bool = True,
                 map_size: int = int(2 ** 40),
                 db_splits: int = 10,
                 index_db_splits: int = 10,
                 **kwargs):
        if host.startswith('~'):
            host = os.path.expanduser(host)
        try:
            os.mkdir(host)
        except FileExistsError:
            pass
        self.host = host
        self.map_size = map_size
        self.db_splits = db_splits
        self.index_db_splits = index_db_splits
        super().__init__(database, host, **kwargs)
        self.has_index = index
        self._setup_dbs()

    def _get_subdb(self, item: str, attr: bool = False) -> str:
        splits = self.index_db_splits if attr else self.db_splits
        return str(hashlib.md5(item.encode()).digest()[0] % splits)

    def _setup_dbs(self) -> None:
        try:
            with open(os.path.join(self.host, 'db_info.json'), 'r') as i:
                db_info = json.load(i)

            if any((
                    db_info['indexed'] != self.has_index,
                    db_info['map_size'] != self.map_size,
                    db_info['db_splits'] != self.db_splits,
                    db_info['index_splits'] != self.index_db_splits,
            )):
                import warnings
                warnings.warn(
                    'LMDB SeqDB settings disagree with on-disk format. Switching to match on-disk.',
                    RuntimeWarning
                )
                self.has_index = db_info['indexed']
                self.map_size = db_info['map_size']
                self.db_splits = db_info['db_splits']
                self.index_db_splits = db_info['index_splits']
        except FileNotFoundError:
            pass

        self.db: Dict[str] = {}
        for i in range(self.db_splits):
            self.db[str(i)] = lmdb.open(os.path.join(self.host, str(i) + '.lmdb'),
                                        map_size=self.map_size // self.db_splits,
                                        writemap=True, map_async=True, readahead=False)
        if self.has_index:
            self.index_dbs: Dict[str] = {}
            for index in self.indices:
                for i in range(self.index_db_splits):
                    self.index_dbs[index + str(i)] = \
                        lmdb.open(os.path.join(self.host, index + str(i) + '.lmdb'),
                                  map_size=self.map_size // self.index_db_splits,
                                  writemap=True, map_async=True, readahead=False)
        with open(os.path.join(self.host, 'db_info.json'), 'w') as o:
            json.dump({'indexed': self.has_index,
                       'map_size': self.map_size,
                       'db_splits': self.db_splits,
                       'index_splits': self.index_db_splits}, o)

    def _reset(self) -> None:
        for db in self.db.values():
            db.close()
        del self.db
        if self.has_index:
            for db in self.index_dbs.values():
                db.close()
            del self.index_dbs
        for subdir in os.listdir(self.host):
            filename = os.path.join(self.host, subdir)
            if os.path.isdir(filename):
                shutil.rmtree(filename)
            else:
                os.remove(filename)
        self._setup_dbs()

    def get_item(self, item: str) -> Union[SeqRecord, None]:
        with self.db[self._get_subdb(item)].begin() as txn:
            t = txn.get(item.encode())
        if not t and self.has_index:
            for attr in (a for a in self.ids if a != '_id'):
                subdb = attr + self._get_subdb(item, True)
                with self.index_dbs[subdb].begin() as txn:
                    t = txn.get(item.encode())
                if t:
                    with self.db[self._get_subdb(t.decode())].begin() as txn:
                        t = txn.get(t)
                        break
        if t is None:
            return None
        return self._extract_seqrecord(t)

    def get_iter(self) -> Generator[SeqRecord, None, None]:
        for i in range(self.db_splits):
            with self.db[str(i)].begin() as txn:
                cursor = txn.cursor()
                for entry in cursor.iternext(keys=False):
                    yield self._extract_seqrecord(entry)

    def get_iterkeys(self) -> Generator[str, None, None]:
        for i in range(self.db_splits):
            with self.db[str(i)].begin() as txn:
                cursor = txn.cursor()
                for entry in cursor.iternext(values=False):
                    yield entry.decode()

    def get_keys(self) -> List[str]:
        return list(self.get_iterkeys())

    def length(self) -> int:
        total = 0
        for i in range(self.db_splits):
            with self.db[str(i)].begin() as txn:
                total += txn.stat()['entries']
        return total

    def get_by(self, attr: str, value: str) -> List[SeqRecord]:
        ret = []
        if attr == '_id':
            ret.append(self.get_item(value))
            return ret
        if self.has_index:
            subdb = attr + self._get_subdb(value, True)
            with self.index_dbs[subdb].begin() as txn:
                db = self.index_dbs[subdb].open_db(dupsort=True)
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
        bpid = pid.encode()
        with self.db[self._get_subdb(pid)].begin(write=True) as txn:
            txn.put(bpid, protein['raw_record'])
        if self.has_index:
            for attr in self.indices:
                if attr in protein:
                    if isinstance(protein[attr], list):
                        for idx in protein[attr]:
                            subdb = attr + self._get_subdb(idx, True)
                            with self.index_dbs[subdb].begin(write=True) as txn:
                                db = self.index_dbs[subdb].open_db(dupsort=True)
                                txn.put(idx.encode(), bpid, db=db)
                    elif isinstance(protein[attr], str):
                        subdb = attr + self._get_subdb(protein[attr], True)
                        with self.index_dbs[subdb].begin(write=True) as txn:
                            db = self.index_dbs[subdb].open_db(dupsort=True)
                            txn.put(protein[attr].encode(), bpid, db=db)
                    else:
                        subdb = attr + self._get_subdb(str(protein[attr]), True)
                        with self.index_dbs[subdb].begin(write=True) as txn:
                            db = self.index_dbs[subdb].open_db(dupsort=True)
                            txn.put(str(protein[attr]).encode(), bpid, db=db)

        return True

    def _add_from_handles(self, handles: Iterable[BinaryIO], filter_fn: Callable = None,
                          total: int = None, loud: bool = False, fake: bool = False) -> None:
        raw_protein_records = itertools.chain(*[parse_raw_swiss(handle, filter_fn) for handle in handles])
        for record in tqdm(raw_protein_records, disable=(not loud), total=total, smoothing=0.1):
            if not fake:
                self.add_protein(self.create_protein_func(record))
