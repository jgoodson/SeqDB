try:
    from itertools import izip_longest as zip_longest
except ImportError:
    from itertools import zip_longest
from datetime import datetime
from collections import defaultdict


def grouper(iterable, n):
    """
    Splits any iterable into chunks of size n

    :param iterable: any iterable
    :param n: size of chunks
    :return: iterable of iterables
    """
    args = [iter(iterable)] * n
    return (_ for _ in zip_longest(*args, fillvalue=None) if not _ is None)


def get_date(record):
    """
    Pulls last update date from Uniprot SeqRecord and returns a datetime

    :param record: SeqRecord object from Uniprot
    :return: datetime of the most recent annotation update
    """
    months = {
        'JAN': 1, 'FEB': 2, 'MAR': 3,
        'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9,
        'OCT': 10, 'NOV': 11, 'DEC': 12,
    }
    day, month, year = record.annotations['date_last_annotation_update'].split('-')
    return datetime(int(year), months[month], int(day))


def get_refs(record):
    """
    Extractions all database references from Uniprot protein SeqRecord

    :param record: Uniprot protein SeqRecord
    :return: dictionary containing all database references in SeqRecord
    """
    result = defaultdict(list)
    result['Uni_name'] = [record.name]
    for ref in record.dbxrefs:
        ref = ref.split(':', 1)
        result[ref[0]].append(ref[1])
    return result

