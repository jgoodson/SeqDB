from collections import defaultdict
from datetime import datetime
from io import StringIO as IOFunc

import zstd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _get_date(dateline: str) -> datetime:
    months = {
        'JAN': 1, 'FEB': 2, 'MAR': 3,
        'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9,
        'OCT': 10, 'NOV': 11, 'DEC': 12,
    }
    day, month, year = dateline.split()[1].strip(',').split('-')
    return datetime(int(year), months[month], int(day))


def _create_record_swiss(raw_record: bytes, compressor: zstd.ZstdCompressor) -> dict:
    lines = raw_record.decode()[:500].split('\n')
    return dict(
        _id=lines[1].split()[1].strip(';'),
        raw_record=compressor.compress(raw_record),
    )


def _create_protein_swiss(raw_record: bytes, compressor: zstd.ZstdCompressor) -> dict:
    lines = raw_record.decode().split('\n')
    desc_lines = []
    refs = defaultdict(list)
    genome = []
    taxid = -1
    dateline = ''
    for line in lines:
        s = line[:2]
        if s == 'CC' or s == '  ':
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

    return dict(
        _id=lines[1].split()[1].strip(';'),
        genome=''.join(genome),
        taxid=taxid,
        description=' '.join(desc_lines),
        updated=_get_date(dateline),
        raw_record=compressor.compress(raw_record),
        Uni_name=[lines[0].split()[1]],
        **refs,
    )


def _extract_seqrecord(raw_record: bytes, decompressor: zstd.ZstdDecompressor) -> SeqRecord:
    return SeqIO.read(IOFunc(decompressor.decompress(raw_record).decode()), 'swiss')
