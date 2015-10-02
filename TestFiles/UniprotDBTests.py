from unittest import TestCase
from cStringIO import StringIO
import tempfile

import UniprotDB


class CreateTest(TestCase):

    def test(self):
        UniprotDB.create_index(['TestFiles/test.dat.bgz'], self.database)
        index = UniprotDB.SeqDB(self.database)
        self.assertEqual(index.get('Q92AT0').id, 'Q92AT0')

    def setUp(self):
        self.tempdb = tempfile.NamedTemporaryFile()
        self.database = 'sqlite:///{}'.format(self.tempdb.name)

    def tearDown(self):
        self.tempdb.close()


class SimpleSqlite(TestCase):

    def setUp(self):
        self.database = 'sqlite:///TestFiles/test.sqlite'
        self.index = UniprotDB.SeqDB(self.database)


class IterationTest(SimpleSqlite):

    def test(self):
        self.assertEqual(self.index.iterkeys().next(), 'Q92AT0')


class KeysTest(SimpleSqlite):

    def test(self):
        self.assertEqual(self.index.keys()[0], 'Q92AT0')


class LenTest(SimpleSqlite):

    def test(self):
        self.assertEqual(len(self.index), 1)
