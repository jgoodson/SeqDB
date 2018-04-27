import unittest
import tempfile

from UniprotDB import UniprotDB

class CreateTest(unittest.TestCase):

    def test_get(self):
        self.assertEqual(self.db.get('Q92AT0').id, 'Q92AT0')

    def test_iter(self):
        self.assertIsInstance(next(self.db.iterkeys()), str)

    def test_keys(self):
        self.assertIn('Q92AT0', self.db.keys())

    def test_len(self):
        self.assertEqual(len(self.db), 900)

    def test_getby(self):
        self.assertEqual(self.db.get_by('_id', 'Q92AT0')[0].id, "Q92AT0")

    def test_fetch(self):
        self.assertEqual(self.db.get('P0AFG0').id, 'P0AFG0')

    def test_update(self):
        with open('TestFiles/test.dat.bgz', 'rb') as h:
            self.db.update([h])

    def setUp(self):
        self.tempdb = tempfile.NamedTemporaryFile()
        self.database = 'test_uni'
        self.db = UniprotDB.create_index(['TestFiles/testbig.dat.gz'], database=self.database)

    def tearDown(self):
        self.tempdb.close()

