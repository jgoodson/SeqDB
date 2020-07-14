import collections
from typing import Callable, Iterable, Union, Generator, List

try:
    from cStringIO import StringIO as IOFunc
except ImportError:
    from io import BytesIO as IOFunc

from Bio.SeqRecord import SeqRecord
import requests
from requests.exceptions import SSLError, ConnectionError

import gzip

sprot_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
trembl_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz'
trembl_taxa_prefix = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions' \
                     '/uniprot_trembl_{}.dat.gz '

query_req = 'https://www.uniprot.org/uniprot/?query={}&format=list'
fetch_req = 'https://www.uniprot.org/uniprot/{}.txt'
uniparc_s_req = 'http://www.uniprot.org/uniparc/?query={}&format=list'
uniparc_f_req = 'http://www.uniprot.org/uniparc/{}.xml'


def search_uniprot(value: str, retries: int = 3) -> Generator[bytes, None, None]:
    possible_ids = []
    for x in range(retries):
        try:
            possible_ids = requests.get(query_req.format(value)).content.split()
            break
        except (SSLError, ConnectionError):
            pass

    raw_record = None
    for pid in possible_ids[:5]:
        for x in range(retries):
            try:
                raw_record = requests.get(fetch_req.format(pid.decode())).content
                break
            except (SSLError, ConnectionError):
                pass
        if raw_record:
            yield raw_record


class SeqDB(collections.abc.Mapping):

    def __init__(self, database: str = 'uniprot',
                 host: Union[str, tuple] = (),
                 dbtype: str = 'lmdb',
                 on_demand: bool = False, **kwargs):
        if dbtype == 'mongo':
            from UniprotDB.MongoDB import MongoDatabase as dbtype
        elif dbtype == 'mongoasync':
            from UniprotDB.AsyncMongoDB import MongoDatabase as dbtype
        elif dbtype == 'lmdb':
            from UniprotDB.LMDB import RawLMDBDatabase as dbtype
        else:
            raise ValueError(f'dbtype: {dbtype} not known')
        self.db = dbtype(database, host, **kwargs)
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


def create_index(flatfiles: Iterable, host: Union[str, tuple] = (), database: str = 'uniprot',
                 filter_fn: Callable[[bytes], bool] = None, **kwargs) -> SeqDB:
    """
    Given a list of SwissProt flatfile filenames in bgz format and a SQLAlchemy database
    identifier, fill the database with the protein entries and returns a SeqDB object.
    """
    s = SeqDB(database, host, **kwargs)
    handles = [gzip.open(f, 'rb') for f in flatfiles]
    s.db.initialize(handles, filter_fn=filter_fn)
    for f in handles:
        f.close()
    return s
