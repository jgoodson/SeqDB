__author__ = 'jonathan'

from sys import stderr
class JEntrez ():

    def __init__(self, email, repo='.'):
        from Bio import Entrez, SeqIO
        from time import sleep,time

        self.Entrez = Entrez
        self.SeqIO = SeqIO
        self.Entrez.email = email
        self.repo = repo
        self.genomes=set()
        self.gmatches=set()
        self.time = time
        self.sleep = sleep
        self.last_ncbi_search = time()

        self.getCache()

    def getCache(self):
        from os import listdir
        files = listdir(self.repo)
        for f in files:
            name = f.rpartition('.')[0]
            match = f.split('.')[0]
            if not name in self.genomes:
                self.genomes.add(name)
                self.gmatches.add(match)

    def getRepoFasta(self, accession_number):
        if accession_number in self.genomes:
            return self.repo+'/'+accession_number+'.fa'
        for x in range(10, 0, -1):
            if accession_number+'.'+str(x) in self.genomes:
                return self.repo+'/'+accession_number+'.'+str(x)+'.fa'

    def getRepoVersionAcc(self, accession_number):
        if accession_number in self.genomes:
            return accession_number
        for x in range(10, 0, -1):
            if '.'.join([accession_number, str(x)]) in self.genomes:
                return '.'.join([accession_number, str(x)])

    def getGenome(self, acc):

        if acc.split('.')[0] in self.gmatches:
            return self.SeqIO.read(self.repo+'/'+self.getRepoVersionAcc(acc)+'.gb', 'gb')

        else:
            if self.time()-self.last_ncbi_search < 1:
                self.sleep(self.time()-self.last_ncbi_search)

        #stderr.write("fetching from NCBI\n")
        net_handle = self.Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
        record = self.SeqIO.read(net_handle, 'gb')
        net_handle.close()

        #self.SeqIO.write(record, ''.join([self.repo,'/',record.id,'.fa']),'fasta')
        self.SeqIO.write(record, ''.join([self.repo,'/',record.id,'.gb']),'gb')

        self.getCache()
        return record

    def get_fasta(gene, gb):
        out=''
        for f in gb.features:
            if 'gene' in f.qualifiers.keys() and f.type=='CDS' and gene in f.qualifiers['gene']:
                s= gb.seq[f.location.start:f.location.end]
                out += '>'+f.qualifiers['gene'][0] + '\n'
                out += str(s) if f.strand == 1 else str(s.reverse_complement())
        return out