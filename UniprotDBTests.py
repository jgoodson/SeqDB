import gzip
import os
import unittest

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

from UniprotDB import UniprotDB
from UniprotDB.SwissProtUtils import filter_proks

ondemand = bool(os.environ.get('TEST_INTERNET'))


class SeqDBTest(object):
    def test_get(self):
        self.assertEqual(self.db.get('Q92AT0').id, 'Q92AT0')

    @unittest.skipUnless(ondemand, "on_demand not enabled for testing")
    def test_get_missing(self):
        self.db.on_demand = True
        self.assertEqual(self.db.get('P0A784').id, 'P0A784')
        self.assertEqual(self.db.get('ORN_HUMAN').id, 'Q9Y3B8')

    # Disabled because it isn't reasonable to hook the get_by() function to download
    # from Uniprot because of the potential huge amount of data it would need to pull
    # @unittest.skipUnless(ondemand, "on_demand not enabled for testing")
    # def test_getby_missing(self):
    #    self.db.on_demand = True
    #    self.assertEqual(self.db.get_by('_id', 'P0A784')[0].id, 'P0A784')

    def test_iter(self):
        self.assertIsInstance(next(self.db.iterkeys()), str)

    def test_iter2(self):
        self.assertIsInstance(next(iter(self.db)).id, str)

    def test_keys(self):
        self.assertIn('Q92AT0', self.db.keys())

    def test_len(self):
        self.assertEqual(len(self.db), 1)

    def test_getby(self):
        self.assertEqual(self.db.get_by('RefSeq', 'WP_010990982.1')[0].id, "Q92AT0")
        self.assertEqual(self.db.get_by('_id', 'Q92AT0')[0].id, "Q92AT0")
        self.assertEqual(self.db.get_by('Uni_name', '12OLP_LISIN')[0].id, "Q92AT0")

    def test_fetch(self):
        self.assertEqual(self.db.get('Q92AT0').id, 'Q92AT0')

    def test_update(self):
        with gzip.open('TestFiles/testbig.dat.gz', 'rb') as h:
            self.db.update([h])
        with gzip.open('TestFiles/testbig.dat.gz', 'rb') as h:
            ids = {line.split()[1].decode() for line in h if line.startswith(b'ID')}
        inserted_ids = {e.name for e in self.db}
        self.assertEqual(inserted_ids, ids)

    def test_update_filtered(self):
        with gzip.open('TestFiles/testbig.dat.gz', 'rb') as h:
            self.db.update([h], filter_fn=filter_proks)
        self.assertEqual(len({e.name for e in self.db}), 70)


@unittest.skipUnless(HAS_MONGO, "requires pymongo")
class MongoTest(unittest.TestCase, SeqDBTest):

    def setUp(self):
        self.database = 'test_uni2'

        import os
        if db_host := os.environ.get('TEST_DB_HOST'):
            self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'],
                                             host=(db_host,),
                                             database=self.database,
                                             dbtype='mongo')
        else:
            self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'],
                                             database=self.database,
                                             dbtype='mongo')

    def tearDown(self):
        pass


@unittest.skipUnless(HAS_MOTOR, "requires motor")
class AsyncTest(unittest.TestCase, SeqDBTest):

    def setUp(self):
        self.database = 'test_uni2'

        import os
        if db_host := os.environ.get('TEST_DB_HOST'):
            self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'], host=(db_host,),
                                             database=self.database, dbtype='mongoasync')
        else:
            self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'],
                                             database=self.database, dbtype='mongoasync')


class LMDBTest(unittest.TestCase, SeqDBTest):

    def setUp(self):
        self.database = 'test_uni2'

        self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'], host='seqdb_test',
                                         database=self.database,
                                         dbtype='lmdb', map_size=int(1024 * 1024 * 1024))

    def tearDown(self):
        import shutil
        for env in self.db.db.db.values():
            env.close()
        for env in self.db.db.index_dbs.values():
            env.close()
        shutil.rmtree('seqdb_test')


if __name__ == '__main__':
    unittest.main()
