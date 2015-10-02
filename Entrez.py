import re
import os
from configparser import ConfigParser
import logging

from Bio import Entrez

from JEntrez import JEntrez

logger = logging.getLogger(__name__)

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
        logger.error('Feature not found {}'.format(protein.id))
        return genome, False
    return genome, feature

