from __future__ import print_function, division

import sys
import itertools
from time import time
from datetime import date
from multiprocessing import Pool
from cStringIO import StringIO
import blosc

from Bio import SeqIO
import sqlalchemy
import sqlalchemy.orm
from sqlalchemy import Column, Integer, String, Date, LargeBinary
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

from SwissProtUtils import parse_raw_swiss, filter_proks

Base = declarative_base()


class UniprotProtein(Base):
    __tablename__ = 'proteins'

    id = Column(String, primary_key=True)

    #dbxrefs
    EMBL = Column(String, index=True)
    RefSeq = Column(String, index=True)
    STRING = Column(String, index=True)
    GeneID = Column(String, index=True)
    PIR = Column(String, index=True)
    Unigene = Column(String, index=True)
    Uni_name = Column(String, index=True)

    #annotation
    genome = Column(String, index=True)
    taxid = Column(Integer, index=True)
    description = Column(String)
    updated = Column(Date, index=True)

    #raw_data
    raw_record = Column(LargeBinary)

    attributes = ['id', 'genome', 'EMBL', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'taxid', 'description', 'updated', 'raw_record']
    indexes = ['id', 'genome', 'EMBL', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'taxid', 'updated']

    def __repr__(self):
        data = {
            'id': self.id,
            'genome': self.genome,
            'EMBL': self.EMBL,
            'RefSeq': self.RefSeq,
            'STRING': self.STRING,
            'GeneID': self.GeneID,
            'PIR': self.PIR,
            'Uni_name': self.Uni_name,
            'taxid': self.taxid,
            'description': self.description,
            'updated': self.updated,
        }
        return '{id} - {genome} ({taxid})\n'\
               'Description={description}' \
               'EMBL={EMBL}\n' \
               'RefSeq={RefSeq}\n' \
               'STRING={STRING}\n' \
               'GeneID={GeneID}\n' \
               'PIR={PIR}\n' \
               'Uni_name={Uni_name}'.format(**data)

def grouper(iterable, n):
    args = [iter(iterable)] * n
    return (_ for _ in itertools.izip_longest(*args, fillvalue=None) if not _ is None)

def get_refs(record):
    good = ['RefSeq', 'STRING', 'GeneID', 'PIR', 'Unigene', 'EMBL']
    result = {'Uni_name': record.name}
    for ref in record.dbxrefs:
        ref = ref.split(':')
        if ref[0] in good:
            result[ref[0]] = ref[1]
    return result


def add_proteins(sequences, database):
    cs = 1000

    p = Pool(15, init_worker, [database])

    chunks = grouper(sequences, cs)

    completed = 0
    current = 0
    t = time()
    start = t
    for done in p.imap_unordered(add_chunk, chunks, 10):
        completed += done
        current += done
        if current >= 50000:
            print('-Completed {} records at {:.2f}r/s ({} total) in {:.1f}s'.format(
                current,
                current/(time()-t),
                completed,
                time()-start), file=sys.stderr)
            current = 0
            t = time()

    p.close()
    p.join()


def _get_date(record):
    months = {
        'JAN': 1, 'FEB': 2, 'MAR': 3,
        'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9,
        'OCT': 10, 'NOV': 11, 'DEC': 12,
    }
    day, month, year = record.annotations['date_last_annotation_update'].split('-')
    return date(int(year), months[month], int(day))


def create_protein(record, raw_record):

    return UniprotProtein(
        id=record.id,
        genome=record.annotations['organism'],
        taxid=record.annotations['ncbi_taxid'][0],
        description=record.description,
        updated=_get_date(record),
        raw_record=blosc.compress(raw_record, 1, cname='zlib', clevel=9),
        **get_refs(record)
    )


def add_chunk(chunk):
    global Session
    session = Session()
    x = 0
    for raw_record in (k for k in chunk if k):
        record = SeqIO.read(StringIO(raw_record), 'swiss')
        protein = create_protein(record, raw_record)
        session.merge(protein)
        x += 1

    session.commit()
    session.close()
    return x

def init_worker(db):
    global Session
    p_engine = sqlalchemy.create_engine(db)
    Session = sqlalchemy.orm.sessionmaker(bind=p_engine)


def initialize_database(seq_files, database, filter_fn=None):
    print("--initializating database\n", file=sys.stderr)
    engine = sqlalchemy.create_engine(database)
    Base.metadata.create_all(engine)

    proteins = itertools.chain(*[parse_raw_swiss(f, filter_fn) for f in seq_files])
    add_proteins(proteins, database)
    print("--initialized database\n", file=sys.stderr)
