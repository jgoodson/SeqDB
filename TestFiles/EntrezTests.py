from unittest import TestCase
from cStringIO import StringIO
import gzip
import tempfile
import shutil

from Bio import SeqIO

import Entrez


class BaseRepo(TestCase):

    def setUp(self):
        self.repo = tempfile.mkdtemp()
        self.e = Entrez.JEntrez('fakeemail@sorryforthetest.edu', self.repo)

    def tearDown(self):
        shutil.rmtree(self.repo)


class FetchTest(BaseRepo):

    def runTest(self):
        genome = self.e.getGenome('NC_000964.3')
        self.assertEqual(genome.id, 'NC_000964.3')


class RepoTest(BaseRepo):

    def runTest(self):
        with_ver = self.e.getGenome('NC_000964.3')
        without_ver = self.e.getGenome('NC_000964')
        self.assertEqual(with_ver.id, without_ver.id)


class GetSourceTest(TestCase):

    def runTest(self):
        p = SeqIO.read(gzip.open('test.dat.bgz'), 'swiss')
        embl = [xref.split(':')[1] for xref in p.dbxrefs if 'EMBL' in xref][-1]
        source_seq, feature = Entrez.get_source_seq(p)
        self.assertIn(embl, source_seq.id)