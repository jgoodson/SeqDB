import tempfile
from unittest import TestCase

from Bio import SeqIO

import Database


class SimpleSqlite(TestCase):

    def setUp(self):
        self.tempdb = tempfile.NamedTemporaryFile()
        self.database = 'sqlite:///{}'.format(self.tempdb.name)
        Database.initialize_database(['test.dat.bgz'], self.database)
        self.SessionMaker = Database.get_sessionmaker(self.database)
        with open('test2.dat') as test_data:
            self.raw_record = test_data.read()

    def tearDown(self):
        pass


class InitTest(SimpleSqlite):

    def runTest(self):
        # Since the setUp() includes initialization, do nothing
        pass


class AddTest(SimpleSqlite):

    def setUp(self):
        super(AddTest, self).setUp()
        self.record = SeqIO.read('test2.dat', 'swiss')
        self.session = self.SessionMaker()

    def tearDown(self):
        self.session.close()
        super(AddTest, self).tearDown()


class AddSingleTest(AddTest):

    def runTest(self):
        Database.add_protein(self.raw_record, self.session)
        added_id = self.session.query(Database.UniprotProtein.id).first()[0]
        self.assertEqual(self.record.id, added_id)


class AddMultTest(AddTest):

    def runTest(self):
        Database.add_proteins_multi([self.raw_record, self.raw_record], self.database, processes=1)
        added_id = self.session.query(Database.UniprotProtein.id).first()[0]
        self.assertEqual(self.record.id, added_id)

class AddMultiprocessingTest(AddTest):

    def runTest(self):
        Database.add_proteins_multi([self.raw_record, self.raw_record], self.database)