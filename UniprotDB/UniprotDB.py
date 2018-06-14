import collections
try:
    from cStringIO import StringIO as IOFunc
except ImportError:
    from io import BytesIO as IOFunc

from UniprotDB.MongoDB import MongoDatabase
import requests
from requests.exceptions import SSLError

query_req = 'https://www.uniprot.org/uniprot/?query={}&format=list'
fetch_req = 'https://www.uniprot.org/uniprot/{}.txt'
uniparc_s_req = 'http://www.uniprot.org/uniparc/?query={}&format=list'
unipart_f_req = 'http://www.uniprot.org/uniparc/{}.xml'

def search_uniprot(value, retries=3):
    for x in range(retries):
        try:
            possible_ids = requests.get(query_req.format(value)).content.split()
            break
        except SSLError:
            pass

    for id in possible_ids[:5]:
        for x in range(retries):
            try:
                raw_record = requests.get(fetch_req.format(id.decode())).content
                break
            except SSLError:
                pass
        if raw_record:
            yield raw_record


class SeqDB(collections.Mapping):


    def __init__(self, database='uniprot', host=(), dbtype=MongoDatabase):
        self.db = dbtype(database, host)
        self.database = database

    def initialize(self, flatfiles, filter):
        self.db.initialize(flatfiles, self.database, filter)

    def __getitem__(self, item):
        r = self.db.get_item(item)
        if not r:
            raw_records = search_uniprot(item)
            for raw_record in raw_records:
                if self.db.add_protein(raw_record, test=item):
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
                if self.db.add_protein(raw_record, test=value, test_attr=attr):
                    r = self.db.get_item(value)
                    break
        return r

    def update(self, handles, filter_fn=None, n_seqs=None, loud=False):
        """
        Update from a Uniprot release file
        :param handle: Streaming handle for a SwissProt format flatfile (gzipped)
        :return: None
        """
        self.db.update(handles, filter_fn=filter_fn, total=n_seqs, loud=loud)



def create_index(flatfiles, host=(), database='uniprot', filter=None):
    """
    Given a list of SwissProt flatfile filenames in bgz format and a SQLAlchemy database
    identifier, fill the database with the protein entries and returns a SeqDB object.
    """
    s = SeqDB(database, host)
    handles = [open(f, 'rb') for f in flatfiles]
    s.db.initialize(handles, filter_fn=filter)
    for f in handles:
        f.close()
    return s


if __name__ == '__main__':
    s = SeqDB()
    print(s['A0A1Q5DVX7'])