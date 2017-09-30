from __future__ import print_function

import collections



class SeqDB(collections.Mapping):

    ids = ['_id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name']

    def __init__(self, database=None):
        self.database = database

    def __getitem__(self, item):
        return self.database.get(item)

    def __iter__(self):
        for item in self.database.iter():
            yield item

    def keys(self):
        return self.database.keys()

    def __len__(self):
        return self.database.len()

    def get_by(self, attr, value):
        return self.database.get_by(attr, value)
