from __future__ import print_function, division

import itertools
import sys
import os.path
import bson

try:
    #Python2
    from itertools import izip_longest as zip_longest
    from itertools import imap
except ImportError:
    #Python3
    from itertools import zip_longest
    imap = map

from multiprocessing import Pool, cpu_count
from io import BytesIO as IOFunc
import zstd

from Bio import SeqIO

import pymongo

from UniprotDB.SwissProtUtils import parse_raw_swiss
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None
from UniprotDB.Utilities import grouper, get_date, get_refs


class MongoDatabase(object):
    ids = ['_id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'EMBL']

    def __init__(self, host, database):
        self.host = host
        self.database = database
        self.client = pymongo.MongoClient(*host)
        self.col = self.client[database].proteins
        self.initialized = False
        try:
            dict_data = self.client[database].proteins.find_one('dict_data')['dict_data']
            self.decompressor = zstd.ZstdDecompressor(dict_data=zstd.ZstdCompressionDict(dict_data))
            self.initialized = True
        except Exception as e:
            print("Database not initialized\n", file=sys.stderr)
        self.decomp = lambda data: self.decompressor.decompress(data)


    def get(self, item):
        """
        Pull a protein record from database given any of the ID types

        :param item: protein ID
        :return: SeqRecord
        """
        t = self.col.find_one({'$or': [{i: item} for i in self.ids]})
        if t is None:
            return None

        r = SeqIO.read(IOFunc(self.decomp(t['raw_record'])), 'swiss')
        return r


    def iter(self):
        for entry in self.col.find({'Uni_name': {'$exists': True}}):
            yield SeqIO.read(IOFunc(self.decomp(entry['raw_record'])), 'swiss')

    def iterkeys(self):
        keys = (i['_id'] for i in self.col.find({'Uni_name': {'$exists': True}},
                                                {'_id': 1}))
        return keys

    def keys(self):
        return self.col.find({'Uni_name': {'$exists': True}},
                             {'_id': 1}).distinct('_id')

    def len(self):
        return self.col.count({'Uni_name': {'$exists': True}})

    def get_by(self, attr, value):
        """
        Fetches protein records with a particular value for a given attribute

        :param attr: Attribute to query
        :param value: Value for attribute
        :return: list of SeqRecords
        """
        res = self.col.find({attr: value}, {'raw_record': 1})
        ret = [SeqIO.read(IOFunc(self.decomp(i['raw_record'])), 'swiss') for i in res]
        return ret

    def initialize_database(self, seq_files,
                            filter_fn=None,
                            zstd_dict_file=os.path.realpath(os.path.dirname(__file__))+'/cc16.zstd_dict',
                            chunksize=500,
                            processes=cpu_count()-1):

        #print("--initializating database\n", file=sys.stderr)
        self.client[self.database].proteins.drop()
        with open(zstd_dict_file, 'rb') as i:
            d = i.read()
            self.client[self.database].proteins.insert_one({'_id': 'dict_data', 'dict_data': bson.binary.Binary(d)})
            dict_data = self.client[self.database].proteins.find_one('dict_data')['dict_data']
            self.decompressor = zstd.ZstdDecompressor(dict_data=zstd.ZstdCompressionDict(dict_data))

        proteins = itertools.chain(*[parse_raw_swiss(f, filter_fn) for f in seq_files])
        _add_proteins(proteins, self.host, self.database, chunksize=chunksize, processes=processes)

        self.client[self.database].proteins.create_index([('genome', 1)])
        indices = ['RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'PDB', 'PDBsum', 'EMBL', 'KEGG',
                   'GO', 'InterPro', 'Pfam', 'TIGRFAMs', 'PROSITE']
        for field in indices:
            self.client[self.database].proteins.create_index(keys=[(field, pymongo.ASCENDING)],
                                                   partialFilterExpression={field: {'$exists': True}})

        #print("--initialized database\n", file=sys.stderr)


def _add_proteins(sequences, host, database, chunksize, processes):
    """
    Function to bulk add protein sequences to MongoDB database

    :param sequences: iterable of SeqRecords
    :param host: MongoDB server hostname
    :param database: MongoDB database
    :param chunksize: number of records to insert per database call
    :param processes: Number of workers
    :return: None
    """

    #p = Pool(processes, _init_worker, (host, database))
    _init_worker(host, database)
    chunks = grouper(sequences, chunksize)

    if tqdm: pbar=tqdm()
    #for done in p.imap_unordered(_add_chunk, chunks):
    for done in imap(_add_chunk, chunks):
        if tqdm: pbar.update(done)
    if tqdm: p.close()


def _init_worker(host=(), database='uniprot'):
    """
    Creates a database client object and compressor object in the current process

    :param host: the MongoDB host information iterable
    :param database: the MongoDB database name
    :return:
    """
    global db, compressor
    client = pymongo.MongoClient(*host)
    db = client[database]

    #prepare Zstd compressor from stored info in database
    dict_data = client[database].proteins.find_one('dict_data')['dict_data']
    compressor = zstd.ZstdCompressor(dict_data=zstd.ZstdCompressionDict(dict_data),
                                     write_content_size=True)


def _add_chunk(chunk):
    """
    Adds an interable of protein record dictionaries created by create_protein() to the currently configured databse.

    :param chunk: Iterable of dicts
    :return: int of number of proteins added
    """

    global db
    collection = db.proteins
    x = 0
    proteins = [_create_protein(raw_record)
                for raw_record in (k for k in chunk if k)]
    x += len(proteins)

    collection.insert_many(proteins)

    return x


def _create_protein(raw_record):
    """
    This function extracts a SwissProt formatted string containing a single complete protein record into a dictionary
    suitable for direct insertion into the database

    :param raw_record: String containing a complete SwissProt format protein entry
    :return: dict containing information about protein record for database entry
    """
    global compressor
    record = SeqIO.read(IOFunc(raw_record), 'swiss')
    return dict(
        _id=record.id,
        genome=record.annotations['organism'],
        taxid=record.annotations['ncbi_taxid'][0],
        description=record.description,
        updated=get_date(record),
        raw_record=bson.binary.Binary(compressor.compress(raw_record)),
        **get_refs(record)
    )