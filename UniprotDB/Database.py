from __future__ import print_function, division

import itertools
import sys
import os.path
from io import BytesIO as IOFunc
from multiprocessing import Pool, cpu_count
from time import time, sleep
import gzip

import sqlalchemy
import sqlalchemy.orm
from Bio import SeqIO
#from progressbar import ProgressBar, Percentage, Bar, AdaptiveETA, RotatingMarker, AdaptiveTransferSpeed
from sqlalchemy import Column, Integer, String, Date, LargeBinary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
import zstd

from UniprotDB.SwissProtUtils import parse_raw_swiss
from UniprotDB.Utilities import grouper, get_date, get_refs


Base = declarative_base()


_known_refs = set(['RefSeq', 'STRING', 'GeneID', 'PIR', 'EMBL'])

class UniprotProtein(Base):
    __tablename__ = 'proteins'

    id = Column(String(25), primary_key=True)

    #dbxrefs
    EMBL = Column(String(20), index=True)
    RefSeq = Column(String(20), index=True)
    STRING = Column(String(50), index=True)
    GeneID = Column(String(20), index=True)
    PIR = Column(String(10), index=True)
    Uni_name = Column(String(20))

    #annotation
    genome = Column(String(255), index=True)
    taxid = Column(Integer, index=True)
    description = Column(String(32768))
    updated = Column(Date, index=True)

    #raw_data
    raw_record = Column(LargeBinary())

    attributes = ['id', 'genome', 'EMBL', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'taxid', 'description', 'updated', 'raw_record']
    indexes = ['id', 'genome', 'EMBL', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'taxid', 'updated']


UP = UniprotProtein

class SQLDatabase(object):
    ids = ['_id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name']

    def __init__(self, database_url):
        self.database = database_url
        p_engine = sqlalchemy.create_engine(self.database)
        self.Session = sqlalchemy.orm.sessionmaker(bind=p_engine)


    def initialize_database(self, seq_files,
                            filter_fn=None,
                            zstd_dict_file=os.path.realpath(os.path.dirname(__file__))+'/cc16.zstd_dict',
                            chunksize=500,
                            processes=cpu_count()-1):

        #print("--initializating database\n", file=sys.stderr)
        engine = sqlalchemy.create_engine(self.database)
        Base.metadata.create_all(engine)

        #TODO find a better place to put the zstd data in the database
        with open(zstd_dict_file, 'rb') as i:
            dict_data = i.read()
        zstd_data = UniprotProtein(
            id='zstd_dict',
            genome=None,
            taxid=None,
            description=None,
            updated=None,
            raw_record=dict_data
        )
        self.decomp = zstd.ZstdDecompressor(dict_data=zstd.ZstdCompressionDict(dict_data))
        session = self.Session()
        session.add(zstd_data)
        session.commit()
        session.close()

        proteins = itertools.chain(*[parse_raw_swiss(f, filter_fn) for f in seq_files])

        _add_proteins(proteins, self.database, chunksize=chunksize, processes=processes)
        #print("--initialized database\n", file=sys.stderr)

    def get(self, item):
        """
        Pull a protein record from database given any of the ID types

        :param item: protein ID
        :return: SeqRecord
        """
        session = self.Session()

        t = session.query(UP.raw_record).filter(
            (UP.id == item) | (UP.STRING == item) | (UP.RefSeq == item) |
            (UP.Uni_name == item) | (UP.EMBL == item) | (UP.GeneID == item)
        ).first()
        session.close()
        if t is None:
            return None
        r = SeqIO.read(IOFunc(self.decomp.decompress(t[0])), 'swiss')

        return r


    def iter(self):
        session = self.Session()
        for row in session.query(UP.raw_record).filter(UP.id!='zstd_dict'):
            yield SeqIO.read(IOFunc(self.decomp.decompress(row[0])), 'swiss')
        session.close()

    def iterkeys(self):
        session = self.Session()
        keys = (i[0] for i in session.query(UP.id).filter(UP.id!='zstd_dict'))
        return keys

    def keys(self):
        session = self.Session()
        keys = [i[0] for i in session.query(UP.id).filter(UP.id!='zstd_dict')]
        session.close()
        return keys

    def len(self):
        session = self.Session()
        size = session.query(UP.Uni_name).filter(UP.id!='zstd_dict').count()
        session.close()
        return size

    def get_by(self, attr, value):
        """
        Fetches protein records with a particular value for a given attribute

        :param attr: Attribute to query
        :param value: Value for attribute
        :return: list of SeqRecords
        """

        session = self.Session()
        res = session.query(UP.raw_record).filter_by(**{attr: value})
        ret = [SeqIO.read(IOFunc(self.decomp.decompress(i[0])), 'swiss') for i in res]
        session.close()
        return ret

def _add_proteins(sequences, database, chunksize, processes):

    p = Pool(processes, _init_worker, [database])
    _init_worker(database)
    chunks = grouper(sequences, chunksize)

    completed = 0

    last_update = time()
    #for done in p.imap_unordered(_add_chunk, chunks):
    for done in map(_add_chunk, chunks):
        completed += done
        if time()-last_update > .25:
            #update pbar
            last_update = time()
        sleep(0.1)
    p.close()
    p.join()


def _add_chunk(chunk):
    global Session
    session = Session()
    x = 0
    proteins = [_create_protein(SeqIO.read(IOFunc(raw_record), 'swiss'), raw_record)
                for raw_record in (k for k in chunk if k)]
    x += len(proteins)
    session.bulk_save_objects(proteins)
    session.commit()
    session.close()

    return x

def _init_worker(db):
    global Session, compressor
    p_engine = sqlalchemy.create_engine(db)
    Session = sqlalchemy.orm.sessionmaker(bind=p_engine)
    session = Session()
    dict_data = session.query(UP.raw_record).filter((UP.id == 'zstd_dict')).first()[0]
    compressor = zstd.ZstdCompressor(dict_data=zstd.ZstdCompressionDict(dict_data),
                                     write_content_size=True)
    session.close()

def _create_protein(record, raw_record):
    global compressor

    return UniprotProtein(
        id=record.id,
        genome=record.annotations['organism'],
        taxid=record.annotations['ncbi_taxid'][0],
        description=record.description,
        updated=get_date(record),
        raw_record=compressor.compress(raw_record),
        **{k:v[0] for k,v in get_refs(record).items() if k in _known_refs}
    )

if __name__ == '__main__':
    d = SQLDatabase('sqlite://'); d.initialize_database('/Users/jonathan/short.dat.gz')