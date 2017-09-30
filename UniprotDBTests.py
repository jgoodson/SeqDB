from __future__ import print_function

import tempfile
import unittest

import pymongo

from UniprotDB.Database import SQLDatabase
from UniprotDB.MongoDB import MongoDatabase
from UniprotDB.SeqDB import SeqDB


class SeqDBTests:
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
    def setUpClass(cls):
        cls.tempdb = tempfile.NamedTemporaryFile()
        cls.sdb = SQLDatabase('sqlite:///' + cls.tempdb.name)
        cls.sdb.initialize_database(['TestFiles/short.dat.gz'])
        cls.s = SeqDB(cls.sdb)

    @classmethod
    def tearDownClass(cls):
        cls.tempdb.close()


class MongoSetup(unittest.TestCase, SeqDBTests):
    def setUp(self):
        if not self.success:
            self.skipTest(reason='MongoDB not available for testing')

    @classmethod
    def setUpClass(cls):
        try:
            pymongo.MongoClient(serverSelectionTimeoutMS=500).server_info()
            cls.success = True
        except pymongo.errors.ServerSelectionTimeoutError:
            cls.success = False

        if cls.success:
            cls.mdb = MongoDatabase((), 'test_uni')
            cls.mdb = MongoDatabase((), 'test_uni')
            cls.mdb.initialize_database(['TestFiles/short.dat.gz'])
            cls.s = SeqDB(cls.mdb)

    @classmethod
    def tearDownClass(cls):
        if cls.success:
            cls.mdb.col.drop()


class SqliteInitTests(unittest.TestCase):
    def setUp(self):
        self.tempdb = tempfile.NamedTemporaryFile()
        self.sdb = SQLDatabase('sqlite:///' + self.tempdb.name)

    def test_single_init(self):
        self.sdb.initialize_database(['TestFiles/short.dat.gz'], chunksize=1, processes=1)

    def multi_init(self):
        self.sdb.initialize_database(['TestFiles/short.dat.gz'], chunksize=1, processes=2)

    def tearDown(self):
        self.tempdb.close()


class MongoInitTests(unittest.TestCase):
    def setUp(self):
        try:
            pymongo.MongoClient(serverSelectionTimeoutMS=500).server_info()
            self.success = True
        except pymongo.errors.ServerSelectionTimeoutError:
            self.success = False
            self.skipTest("MongoDB not available")

        if self.success:
            self.mdb = MongoDatabase((), 'test_uni')

    def test_single_init(self):
        self.mdb.initialize_database(['TestFiles/short.dat.gz'], chunksize=1, processes=1)

    def multi_init(self):
        self.mdb.initialize_database(['TestFiles/short.dat.gz'], chunksize=1, processes=2)


# TODO Test filter_proks

if __name__ == '__main__':
    unittest.main()
