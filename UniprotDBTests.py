from __future__ import print_function

import unittest
import tempfile
import os
import pymongo

from UniprotDB.SeqDB import SeqDB
from UniprotDB.MongoDB import MongoDatabase
from UniprotDB.Database import SQLDatabase

class SeqDBTests():

    def iter_test(self):
        self.assertEqual(len(next(iter(self.s)).id), 6)

    def keys_test(self):
        self.assertEqual(len(self.s.keys()), 4)

    def len_test(self):
        self.assertEqual(len(self.s), 4)

    def test_get_id(self):
        self.assertEqual(self.s['Q6GZX4'].id, 'Q6GZX4')

    def test_get_embl_id(self):
        self.assertEqual(self.s['AY548484'].id, 'Q6GZX4')

    def test_get_refseq_id(self):
        self.assertEqual(self.s['YP_031579.1'].id, 'Q6GZX4')

    def test_get_geneid_id(self):
        self.assertEqual(self.s['2947773'].id, 'Q6GZX4')

    def test_getby(self):
        self.assertEqual(self.s.get_by('taxid', '654924')[0].id, 'Q6GZX4')


class SqliteSetup(unittest.TestCase, SeqDBTests):

    @classmethod
    def setUpClass(self):
        self.tempdb = tempfile.NamedTemporaryFile()
        sdb = SQLDatabase('sqlite:///'+self.tempdb.name)
        sdb.initialize_database(['TestFiles/short.dat.gz'])
        self.s = SeqDB(sdb)

    @classmethod
    def tearDownClass(self):
        self.tempdb.close()

class MongoSetup(unittest.TestCase, SeqDBTests):

    def setUp(self):
        if not self.success:
            self.skipTest(reason='MongoDB not available for testing')

    @classmethod
    def setUpClass(self):
        try:
            pymongo.MongoClient(serverSelectionTimeoutMS=500).server_info()
            self.success = True
        except pymongo.errors.ServerSelectionTimeoutError:
            self.success = False

        if self.success:
            self.mdb = MongoDatabase((), 'test_uni')
            self.mdb = MongoDatabase((), 'test_uni')
            self.mdb.initialize_database(['TestFiles/short.dat.gz'])
            self.s = SeqDB(self.mdb)

    @classmethod
    def tearDownClass(self):
        if self.success:
            self.mdb.col.drop()


if __name__ == '__main__':
    unittest.main()
