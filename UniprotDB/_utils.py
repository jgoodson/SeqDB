import hashlib
from collections import defaultdict
from datetime import datetime
from io import StringIO as IOFunc
from typing import Generator

import requests
import zstandard
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from requests.exceptions import SSLError, ConnectionError

query_req = 'https://www.uniprot.org/uniprot/?query={}&format=list'
fetch_req = 'https://www.uniprot.org/uniprot/{}.txt'
uniparc_s_req = 'http://www.uniprot.org/uniparc/?query={}&format=list'
uniparc_f_req = 'http://www.uniprot.org/uniparc/{}.xml'


def _get_date(dateline: str) -> datetime:
    months = {
        'JAN': 1, 'FEB': 2, 'MAR': 3,
        'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9,
        'OCT': 10, 'NOV': 11, 'DEC': 12,
    }
    day, month, year = dateline.split()[1].strip(',').split('-')
    return datetime(int(year), months[month], int(day))


def _create_protein_swiss(raw_record: bytes, compressor: zstandard.ZstdCompressor) -> dict:
    lines = raw_record.decode().split('\n')
    desc_lines = []
    refs = defaultdict(list)
    genome = []
    in_seq = False
    seq_lines = []
    taxid = -1
    dateline = ''
    for line in lines:
        s = line[:2]
        if s == 'CC':
            continue
        elif s == 'DT':
            dateline = line
        elif s == 'DE':
            desc_lines.append(line.split(maxsplit=1)[1])
        elif s == 'OS':
            genome.append(line.split(maxsplit=1)[1].strip('. '))
        elif s == 'OX':
            taxid = int(line.split('=')[1].split()[0].strip(';'))
        elif s == 'DR':
            ref = line.split(maxsplit=1)[1]
            dec = ref.strip('.').split(';')
            refs[dec[0]].append(dec[1].strip())
        elif in_seq and s == '  ':
            seq_lines.append(line)
        elif s == 'SQ':
            in_seq = True

    seq = ''.join(''.join(seq_lines).split())

    return dict(
        _id=lines[1].split()[1].strip(';'),
        genome=''.join(genome),
        taxid=taxid,
        description=' '.join(desc_lines),
        updated=_get_date(dateline),
        raw_record=compressor.compress(raw_record),
        seq_sha1=hashlib.sha1(seq.encode()).hexdigest(),
        Uni_name=[lines[0].split()[1]],
        **refs,
    )


def _extract_seqrecord(raw_record: bytes, decompressor: zstandard.ZstdDecompressor) -> SeqRecord:
    return SeqIO.read(IOFunc(decompressor.decompress(raw_record).decode()), 'swiss')


def search_uniprot(value: str, retries: int = 3) -> Generator[bytes, None, None]:
    possible_ids = []
    for x in range(retries):
        try:
            possible_ids = requests.get(query_req.format(value)).content.split()
            break
        except (SSLError, ConnectionError):
            pass

    raw_record = None
    for pid in possible_ids[:5]:
        for x in range(retries):
            try:
                raw_record = requests.get(fetch_req.format(pid.decode())).content
                break
            except (SSLError, ConnectionError):
                pass
        if raw_record:
            yield raw_record
