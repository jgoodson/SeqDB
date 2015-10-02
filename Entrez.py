import re
import os
from configparser import ConfigParser

from Bio import Entrez

from JEntrez import JEntrez

conf = ConfigParser()
base = os.path.dirname(os.path.realpath(__file__))
conf.read(base+'/settings.conf')


Entrez.email = conf.get('Entrez', 'email')
repo_location = conf.get('Entrez', 'repo_location')
e = JEntrez(Entrez.email, repo_location)
rs_acc = re.compile('[A-Z]+\d+\.\d')

def get_source_seq(protein):
    embl = [xref.split(':')[1] for xref in protein.dbxrefs if 'EMBL' in xref][-1]
    genome = e.getGenome(embl)
    protein_seq = str(protein.seq)
    for feature in genome.features:
        if feature.type == 'CDS' and 'translation' in feature.qualifiers and feature.qualifiers['translation'][0] == protein_seq:
            break
    else:
        print 'Feature not found {}'.format(protein.id)
        return False
    return genome, feature

