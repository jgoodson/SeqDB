import collections

try:
    from cStringIO import StringIO as IOFunc
except ImportError:
    from io import BytesIO as IOFunc

from UniprotDB.MongoDB import MongoDatabase
import requests
from requests.exceptions import SSLError, ConnectionError

sprot_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
trembl_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz'
trembl_taxa_prefix = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_{}.dat.gz'

query_req = 'https://www.uniprot.org/uniprot/?query={}&format=list'
fetch_req = 'https://www.uniprot.org/uniprot/{}.txt'
uniparc_s_req = 'http://www.uniprot.org/uniparc/?query={}&format=list'
unipart_f_req = 'http://www.uniprot.org/uniparc/{}.xml'


def search_uniprot(value, retries=3):
    possible_ids = []
    for x in range(retries):
        try:
            possible_ids = requests.get(query_req.format(value)).content.split()
            break
        except (SSLError, ConnectionError) as e:
            pass

    for id in possible_ids[:5]:
        for x in range(retries):
            try:
                raw_record = requests.get(fetch_req.format(id.decode())).content
                break
            except (SSLError, ConnectionError) as e:
                pass
        if raw_record:
            yield raw_record


class SeqDB(collections.Mapping):

    def __init__(self, database='uniprot', host=(), dbtype=MongoDatabase, on_demand=False):
        self.db = dbtype(database, host)
        self.database = database
        self.on_demand = on_demand

    def initialize(self, flatfiles, filter):
        self.db.initialize(flatfiles, self.database, filter)

    def __getitem__(self, item):
        r = self.db.get_item(item)
        if not r and self.on_demand:
            raw_records = search_uniprot(item)
            for raw_record in raw_records:
                if self.db.add_record(raw_record, test=item):
                    r = self.db.get_item(item)
                    break
        return r

    def __iter__(self):
        return self.db.get_iter()

    def iterkeys(self):
        return self.db.get_iterkeys()

    def keys(self):
        return self.db.get_keys()

    def __len__(self):
        return self.db.length()

    def get_by(self, attr, value):
        r = self.db.get_by(attr, value)
        if not r:
            raw_records = search_uniprot(value)
            for raw_record in raw_records:
                if self.db.add_record(raw_record, test=value, test_attr=attr):
                    r = self.db.get_item(value)
                    break
        return r

    def update(self, handles, filter_fn=None, n_seqs=None, loud=False, processes=1):
        self.db.update(handles, filter_fn=filter_fn, total=n_seqs, loud=loud, processes=processes)

    def update_swissprot(self, filter_fn=None, processes=1, loud=True):
        import urllib.request
        sprot = urllib.request.urlopen(sprot_url)
        self.update([sprot], filter_fn=filter_fn, loud=loud, processes=processes)
        sprot.close()

    def update_trembl_taxa(self, taxa, filter_fn=None, processes=1, loud=True):
        import urllib.request
        for taxon in taxa:
            taxon_handle = urllib.request.urlopen(trembl_taxa_prefix.format(taxon))
            print("Updating {}".format(taxon))
            self.update([taxon_handle], filter_fn, loud, processes)
            taxon_handle.close()

    def update_trembl_prok(self, filter_fn=None, processes=1, loud=True):
        self.update_trembl_taxa(['bacteria', 'archaea'],
                                filter_fn, processes, loud)

    def update_trembl_euk(self, filter_fn=None, processes=1, loud=True):
        self.update_trembl_taxa(['fungi', 'human', 'invertebrate', 'mammal', 'plant', 'rodent', 'vertebrate', 'virus'],
                                filter_fn, processes, loud)

    def update_trembl(self, filter_fn=None, processes=1, loud=True):
        import urllib.request
        trembl = urllib.request.urlopen(trembl_url)
        self.update([trembl], filter_fn=filter_fn, loud=loud, processes=processes)
        trembl.close()


def create_index(flatfiles, host=(), database='uniprot', filter=None, **kwargs):
    """
    Given a list of SwissProt flatfile filenames in bgz format and a SQLAlchemy database
    identifier, fill the database with the protein entries and returns a SeqDB object.
    """
    s = SeqDB(database, host, **kwargs)
    handles = [open(f, 'rb') for f in flatfiles]
    s.db.initialize(handles, filter_fn=filter)
    for f in handles:
        f.close()
    return s


if __name__ == '__main__':
    s = SeqDB()
    print(s['A0A1Q5DVX7'])
