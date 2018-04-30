from __future__ import print_function

import collections
try:
    from cStringIO import StringIO as IOFunc
except ImportError:
    from io import BytesIO as IOFunc

from MongoDB import MongoDatabase
import requests

from SwissProtUtils import filter_proks

query_req = 'https://www.uniprot.org/uniprot/?query={}&format=list'
fetch_req = 'https://www.uniprot.org/uniprot/{}.txt'

class SeqDB(collections.Mapping):


    def __init__(self, database='uniprot', host=(), dbtype=MongoDatabase):
        self.db = dbtype(database, host)
        self.database = database

    def initialize(self, flatfiles, filter):
        self.db.initialize(flatfiles, self.database, filter)

    def __getitem__(self, item):
        r = self.db.get_item(item)
        if not r:
            possible_ids = requests.get(query_req.format(item)).content.split()
            for id in possible_ids[:5]:
                raw_record = requests.get(fetch_req.format(id.decode())).content
                if self.db.add_protein(raw_record, test=item):
                    break
            r = self.db.get_item(item)
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
            possible_ids = requests.get(query_req.format(value)).content.split()
            for id in possible_ids[:5]:
                raw_record = requests.get(fetch_req.format(id.decode())).content
                if self.db.add_protein(raw_record, test=value, test_attr=attr):
                    break
            r = self.db.get_item(value)
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
    s.db.initialize((open(f, 'rb') for f in flatfiles), database=database, filter_fn=filter)
    return s


if __name__ == '__main__':
    s = SeqDB()
    print(s['A0A1Q5DVX7'])