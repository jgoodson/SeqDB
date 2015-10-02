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

    def tearDown(self):
        pass

class InitTest(SimpleSqlite):

    def runTest(self):
        # Since the setUp() includes initialization, do nothing
        pass

class AddTest(SimpleSqlite):

    def runTest(self):
        with open('test2.dat') as test_data:
            raw_record = test_data.read()
        session = self.SessionMaker()
        Database.add_protein(raw_record, session)

        record = SeqIO.read('test2.dat', 'swiss')
        added_id = session.query(Database.UniprotProtein.id).first()[0]
        self.assertEqual(record.id, added_id)
        session.close()

class AddMultTest(SimpleSqlite):

    def runTest(self):
        with open('test2.dat') as test_data:
            raw_record = test_data.read()
        Database.add_proteins([raw_record, raw_record], self.database)
        record = SeqIO.read('test2.dat', 'swiss')
        session = self.SessionMaker()
        added_id = session.query(Database.UniprotProtein.id).first()[0]
        self.assertEqual(record.id, added_id)
        session.close()