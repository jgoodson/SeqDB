import tempfile
import unittest

from UniprotDB import UniprotDB
from UniprotDB.SwissProtUtils import filter_proks


class CreateTest(unittest.TestCase):

    def test_direct(self):
        self.assertIn('Q92AT0', self.db.db.loop.run_until_complete(self.db.db.col.distinct('_id')))

    def test_get(self):
        self.assertEqual(self.db.get('Q92AT0').id, 'Q92AT0')

    def test_get_missing(self):
        self.assertEqual(self.db.get('P0A784').id, 'P0A784')
        self.assertEqual(self.db.get('ORN_HUMAN').id, 'Q9Y3B8')

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

    def test_fetch(self):
        self.assertEqual(self.db.get('Q92AT0').id, 'Q92AT0')

    def test_update(self):
        with open('TestFiles/testbig.dat.gz', 'rb') as h:
            self.db.update([h])
        self.assertEqual(len(self.db), 900)
        if len(self.db) != 900:
            print(self.db.keys())

    def test_update_filtered(self):
        with open('TestFiles/testbig.dat.gz', 'rb') as h:
            self.db.update([h], filter_fn=filter_proks)
        self.assertEqual(len(self.db), 69)

    def setUp(self):
        self.tempdb = tempfile.NamedTemporaryFile()
        self.database = 'test_uni'
        self.db = UniprotDB.create_index(['TestFiles/test.dat.bgz'], database=self.database)

    def tearDown(self):
        self.tempdb.close()

if __name__ == '__main__':
    unittest.main()