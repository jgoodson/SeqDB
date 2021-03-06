import collections
from typing import Callable, Iterable, Union, Generator, List

try:
    import motor

    HAS_MOTOR = True
except ImportError:
    HAS_MOTOR = False
try:
    import pymongo

    HAS_MONGO = True
except ImportError:
    HAS_MONGO = False

from UniprotDB._utils import search_uniprot

try:
    from cStringIO import StringIO as IOFunc
except ImportError:
    from io import BytesIO as IOFunc

from Bio.SeqRecord import SeqRecord

import gzip

sprot_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
trembl_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz'
trembl_taxa_prefix = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions' \
                     '/uniprot_trembl_{}.dat.gz '


class SeqDB(collections.abc.Mapping):

    def __init__(self, database: str = 'uniprot',
                 host: Union[str, tuple] = '',
                 dbtype: str = 'lmdb',
                 on_demand: bool = False, **kwargs):
        if dbtype == 'mongo':
            if HAS_MONGO:
                from UniprotDB.MongoDB import MongoDatabase as BaseDB
            else:
                raise ModuleNotFoundError('Missing pymongo')
        elif dbtype == 'mongoasync':
            if HAS_MOTOR:
                from UniprotDB.AsyncMongoDB import MongoDatabase as BaseDB
            else:
                raise ModuleNotFoundError('Missing motor')
        elif dbtype == 'lmdb':
            from UniprotDB.LMDB import RawLMDBDatabase as BaseDB
        else:
            raise ValueError(f'BaseDB: {dbtype} not known')
        if host:
            self.db = BaseDB(database, host, **kwargs)
        else:
            self.db = BaseDB(database, **kwargs)
        self.database = database
        self.on_demand = on_demand

    def initialize(self, flatfiles: Iterable, *args, **kwargs) -> None:
        self.db.initialize(flatfiles, *args, **kwargs)

    def __getitem__(self, item: str) -> SeqRecord:
        r = self.db.get_item(item)
        if not r and self.on_demand:
            raw_records = search_uniprot(item)
            for raw_record in raw_records:
                if self.db.add_record(raw_record, test=item):
                    r = self.db.get_item(item)
                    break
        return r

    def __iter__(self) -> Generator[SeqRecord, None, None]:
        return self.db.get_iter()

    def iterkeys(self) -> Generator[str, None, None]:
        return self.db.get_iterkeys()

    def keys(self) -> List[str]:
        return self.db.get_keys()

    def __len__(self) -> int:
        return self.db.length()

    def get_by(self, attr: str, value: str) -> List[SeqRecord]:
        return self.db.get_by(attr, value)

    def update(self, handles: Iterable, filter_fn: Callable = None,
               n_seqs: int = None, loud: bool = False, workers: int = 1) -> None:
        self.db.update(handles, filter_fn=filter_fn, total=n_seqs, loud=loud, workers=workers)

    def update_swissprot(self, filter_fn: Callable[[bytes], bool] = None, workers: int = 1, loud: bool = True) -> None:
        import urllib.request
        sprot = gzip.open(urllib.request.urlopen(sprot_url))
        self.update([sprot], filter_fn=filter_fn, loud=loud, workers=workers)
        sprot.close()

    def update_trembl_taxa(self, taxa: Iterable, filter_fn: Callable[[bytes], bool] = None,
                           workers: int = 1, loud: bool = True) -> None:
        import urllib.request
        for taxon in taxa:
            taxon_handle = gzip.open(urllib.request.urlopen(trembl_taxa_prefix.format(taxon)))
            print("Updating {}".format(taxon))
            self.update([taxon_handle], filter_fn, loud, workers=workers)
            taxon_handle.close()

    def update_trembl_prok(self, filter_fn: Callable[[bytes], bool] = None,
                           workers: int = 1, loud: bool = True) -> None:
        self.update_trembl_taxa(['bacteria', 'archaea'],
                                filter_fn, workers, loud)

    def update_trembl_euk(self, filter_fn: Callable[[bytes], bool] = None,
                          workers: int = 1, loud: bool = True) -> None:
        self.update_trembl_taxa(['fungi', 'human', 'invertebrate', 'mammal', 'plant', 'rodent', 'vertebrate', 'virus'],
                                filter_fn, workers, loud)

    def update_trembl(self, filter_fn: Callable[[bytes], bool] = None,
                      workers: int = 1, loud: bool = True) -> None:
        import urllib.request
        trembl = gzip.open(urllib.request.urlopen(trembl_url))
        self.update([trembl], filter_fn=filter_fn, loud=loud, workers=workers)
        trembl.close()


def create_index(flatfiles: Iterable, host: Union[str, tuple] = (),
                 dbtype: str = 'lmdb', n_jobs: int = 1, **kwargs) -> SeqDB:
    """
    Given a list of SwissProt flatfile filenames in uncompressed, gzip, or zstd format and a database
    host/filename + dbtype, fill the database with the protein entries and returns a SeqDB object.
    """
    from .data_loader import process_main
    s = process_main(flatfiles, host,
                     dbtype=dbtype, initialize=True, verbose=True, n_jobs=n_jobs, **kwargs)

    return s
