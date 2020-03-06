from datetime import datetime
from collections import defaultdict

from Bio import SeqIO
from io import StringIO as IOFunc

def _get_date(dateline):
    months = {
        'JAN': 1, 'FEB': 2, 'MAR': 3,
        'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9,
        'OCT': 10, 'NOV': 11, 'DEC': 12,
    }
    day, month, year = dateline.split()[1].strip(',').split('-')
    return datetime(int(year), months[month], int(day))


def _create_protein_swiss(raw_record, compressor):
    lines = raw_record.decode().split('\n')
    desc_lines = []
    refs = defaultdict(list)
    genome = []
    for l in lines:
        s = l[:2]
        if s == 'CC' or s == '  ':
            continue
        elif s == 'DT':
            dateline = l
        elif s == 'DE':
            desc_lines.append(l.split(maxsplit=1)[1])
        elif s == 'OS':
            genome.append(l.split(maxsplit=1)[1].strip('. '))
        elif s == 'OX':
            taxid = int(l.split('=')[1].split()[0].strip(';'))
        elif s == 'DR':
            ref = l.split(maxsplit=1)[1]
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


def _extract_seqrecord(raw_record, decompressor):
    return SeqIO.read(IOFunc(decompressor.decompress(raw_record).decode()), 'swiss')
