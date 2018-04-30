from __future__ import print_function, division

import sys
import itertools
try:
    from itertools import izip_longest as zip_longest
except ImportError:
    from itertools import zip_longest
import os
from datetime import datetime
from multiprocessing import Pool
from io import BytesIO as IOFunc
import zstd
from collections import defaultdict

from Bio import SeqIO

import pymongo

from SwissProtUtils import parse_raw_swiss, filter_proks
from tqdm import tqdm


def _grouper(iterable, n):
    args = [iter(iterable)] * n
    return (_ for _ in zip_longest(*args, fillvalue=None) if not _ is None)


def _add_proteins(sequence_records, host, database, loud=False):
    cs = 500
    p = Pool(16, _init_worker, (host, database))
    _init_worker(host, database)

    record_chunks = _grouper(sequence_records, cs)
    if loud:
        with tqdm(total=81000000) as pbar:
            for done in p.imap_unordered(_add_chunk, record_chunks):
                pbar.update(done)
    else:
        for done in p.imap_unordered(_add_chunk, record_chunks):
            pass
    p.close()

def _init_worker(host=(), database='uniprot'):
    global collection, compressor
    client = pymongo.MongoClient(*host)
    collection = client[database]
    dict_data = client[database].proteins.find_one('dict_data')['dict_data']
    compressor = zstd.ZstdCompressor(dict_data=zstd.ZstdCompressionDict(dict_data),
                                     write_content_size=True)

def _get_date(dateline):
    months = {
        'JAN': 1, 'FEB': 2, 'MAR': 3,
        'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9,
        'OCT': 10, 'NOV': 11, 'DEC': 12,
    }
    day, month, year = dateline.split()[1].strip(',').split('-')
    return datetime(int(year), months[month], int(day))

def _create_protein(raw_record, compressor):
    lines = raw_record.decode().split('\n')
    desc_lines = []
    refs = defaultdict(list)
    for l in lines:
        s = l[:2]
        if s == 'CC' or s == '  ':
            continue
        elif s == 'DT':
            dateline = l
        elif s == 'DE':
            desc_lines.append(l.split(maxsplit=1)[1])
        elif s == 'OS':
            genome = l.split(maxsplit=1)[1].strip('. ')
        elif s == 'OX':
            taxid = l.split('=')[1].strip(';')
        elif s == 'DR':
            ref = l.split(maxsplit=1)[1]
            dec = ref.strip('.').split(';')
            refs[dec[0]].append(dec[1].strip())


    return dict(
        _id=lines[1].split()[1].strip(';'),
        genome=genome,
        taxid=taxid,
        description=' '.join(desc_lines),
        updated=_get_date(dateline),
        raw_record=compressor.compress(raw_record),
        Uni_name=[lines[0].split()[1]],
        **refs,
    )

def _add_chunk(record_chunk):
    global client, collection, compressor

    x = 0
    proteins = [_create_protein(raw_record, compressor)
                for raw_record in (k for k in record_chunk if k)]
    x += len(proteins)
    collection.proteins.insert_many(proteins)

    return x


class MongoDatabase(object):

    ids = ['_id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name']

    def __init__(self, database, host):
        self.host = host
        self.client = pymongo.MongoClient(*host)
        self.database = database
        try:
            self.col = self.client[database].proteins
            dict_data = self.client[database].proteins.find_one('dict_data')['dict_data']
            self.compressor = zstd.ZstdCompressor(dict_data=zstd.ZstdCompressionDict(dict_data),
                                                  write_content_size=True)
            self.decomp = zstd.ZstdDecompressor(dict_data=zstd.ZstdCompressionDict(dict_data))
        except:
            pass #Do initialization?


    def get_item(self, item):
        t = self.col.find_one({'$or': [{i: item} for i in self.ids]})
        if t is None:
            return None
        r = SeqIO.read(IOFunc(self.decomp.decompress(t['raw_record'])), 'swiss')
        return r


    def get_iter(self):
        for entry in self.col.find({'Uni_name': {'$exists': True}}):
            yield SeqIO.read(IOFunc(self.decomp.decompress(entry['raw_record'])), 'swiss')


    def get_iterkeys(self):
        keys = (i['_id'] for i in self.col.find({'Uni_name': {'$exists': True}},
                                                {'_id': 1}))
        return keys


    def get_keys(self):
        return self.col.distinct('_id')


    def length(self):
        return self.col.count({'Uni_name': {'$exists': True}})


    def get_by(self, attr, value):
        res = self.col.find({attr: value}, {'raw_record': 1})
        ret = [SeqIO.read(IOFunc(self.decomp.decompress(i['raw_record'])), 'swiss') for i in res]
        return ret


    def initialize(self, seq_handles, database='uniprot', filter_fn=None, loud=False):

        if loud:
            print("--initializating database\n", file=sys.stderr)
        self.client[database].proteins.drop()
        with open(os.path.dirname(__file__)+"/cc16.zstd_dict", 'rb') as i:
            d = i.read()
            self.client[database].proteins.insert_one({'_id': 'dict_data', 'dict_data': d})

        raw_protein_recordss = itertools.chain(*[parse_raw_swiss(handle, filter_fn) for handle in seq_handles])

        _add_proteins(raw_protein_recordss, self.host, database, loud=loud)

        self.client[database].proteins.create_index([('genome', 1)])
        indices = ['RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'PDB', 'PDBsum', 'EMBL', 'KEGG',
                   'GO', 'InterPro', 'Pfam', 'TIGRFAMs', 'PROSITE']
        for field in indices:
            self.client[database].proteins.create_index(keys=[(field, pymongo.ASCENDING)],
                                                   partialFilterExpression={field: {'$exists': True}})
        if loud:
            print("--initialized database\n", file=sys.stderr)

    def add_protein(self, raw_record, update=False, test=None, test_attr=None):
        protein = _create_protein(raw_record, self.compressor)

        if test:
            good = False
            if test == protein['_id']:
                good = True
            for ref in [test_attr] if test_attr else self.ids:
                if test in getattr(protein, ref, []):
                    good = True
            if not good:
                return False

        if update:
            self.col.replace_one({'_id': protein['_id']}, protein, upsert=True)
        else:
            self.col.insert_one(protein)
        return True

    def update(self, handles, filter_fn=None, loud=False, total=None):
        raw_protein_records = itertools.chain(*[parse_raw_swiss(handle, filter_fn, check_date=True) for handle in handles])
        added = 0
        with tqdm(total=total) as pbar:
            for record, date in (r for r in raw_protein_records):
                if date > getattr(self.col.find_one(record.split(maxsplit=2)[1]), 'updated', datetime.min):
                    self.add_protein(record, update=True)
                    added += 1
                pbar.update()

