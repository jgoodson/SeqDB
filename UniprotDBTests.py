import gzip
import unittest

from UniprotDB import UniprotDB
from UniprotDB.SwissProtUtils import filter_proks


class MongoTest(unittest.TestCase):

    def test_get(self):
        self.assertEqual(self.db.get('Q92AT0').id, 'Q92AT0')

    def test_get_missing(self):
        if self.ondemand:
            self.assertEqual(self.db.get('P0A784').id, 'P0A784')
            self.assertEqual(self.db.get('ORN_HUMAN').id, 'Q9Y3B8')
        else:
            self.skipTest('Not testing on-demand feature')

    def test_getby_missing(self):
        if self.ondemand:
            self.assertEqual(self.db.get_by('_id', 'P0A784').id, 'P0A784')
        else:
            self.skipTest('Not testing on-demand feature')

    def test_iter(self):
        self.assertIsInstance(next(self.db.iterkeys()), str)

    def test_iter2(self):
        self.assertIsInstance(next(iter(self.db)).id, str)

    def test_keys(self):
        self.assertIn('Q92AT0', self.db.keys())

    def test_len(self):
        self.assertEqual(len(self.db), 1)

    def test_getby(self):
        self.assertEqual(self.db.get_by('_id', 'Q92AT0')[0].id, "Q92AT0")
        self.assertEqual(self.db.get_by('Uni_name', '12OLP_LISIN')[0].id, "Q92AT0")
        self.assertEqual(self.db.get_by('RefSeq', 'WP_010990982.1')[0].id, "Q92AT0")

    def test_fetch(self):
        self.assertEqual(self.db.get('Q92AT0').id, 'Q92AT0')

    def test_update(self):
        with gzip.open('TestFiles/testbig.dat.gz', 'rb') as h:
            self.db.update([h])
        with gzip.open('TestFiles/testbig.dat.gz', 'rb') as h:
            ids = set(l.split()[1].decode() for l in h if l.startswith(b'ID'))
        inserted_ids = set(e.name for e in self.db)
        self.assertEqual(inserted_ids, ids)

    def test_update_filtered(self):
        with gzip.open('TestFiles/testbig.dat.gz', 'rb') as h:
            self.db.update([h], filter_fn=filter_proks)
        self.assertEqual(len(set(e.name for e in self.db)), 70)

    def setUp(self):
        self.database = 'test_uni2'

        import os
        self.ondemand = os.environ.get('TEST_INTERNET')
        db_host = os.environ.get('TEST_DB_HOST')

        if db_host:
            self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'],
                                             host=(db_host,),
                                             database=self.database,
                                             on_demand=True)
        else:
            self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'],
                                             database=self.database,
                                             on_demand=True)

    def tearDown(self):
        pass


class AsyncTest(MongoTest):

    def setUp(self):
        from UniprotDB.AsyncMongoDB import MongoDatabase
        self.database = 'test_uni2'

        import os
        self.ondemand = os.environ.get('TEST_INTERNET')
        db_host = os.environ.get('TEST_DB_HOST')

        if db_host:
            self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'], host=(db_host,),
                                             database=self.database, on_demand=True, dbtype=MongoDatabase)
        else:
            self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'],
                                             database=self.database, on_demand=True, dbtype=MongoDatabase)

class LMDBTest(MongoTest):

    def setUp(self):
        from UniprotDB.LMDB import LMDBDatabase
        import os
        self.database = 'test_uni2'
        self.ondemand = os.environ.get('TEST_INTERNET')

        self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'], host='test.lmdb',
                                         database=self.database, on_demand=True, dbtype=LMDBDatabase)

    def tearDown(self):
        import shutil
        shutil.rmtree('test.lmdb')

if __name__ == '__main__':
    unittest.main()
