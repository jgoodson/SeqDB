import gzip
import re
from datetime import datetime

def _get_date(line):
    months = {
        'JAN': 1, 'FEB': 2, 'MAR': 3,
        'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9,
        'OCT': 10, 'NOV': 11, 'DEC': 12,
    }

    day, month, year = line.decode().split()[1].strip(',').split('-')
    return datetime(int(year), months[month], int(day))

def _get_record(handle, check_date):
    """
    Returns the next complete SwissProt entry in the input handle
    """
    lines = []
    date = None
    for line in handle:
        lines.append(line)
        if check_date and line.startswith(b'DT'):
            date = _get_date(line)
        if line.startswith(b'//'):
            if check_date:
                yield b''.join(lines), date
            else:
                yield b''.join(lines)
            lines = []
            date = None


def filter_proks(record):
    """
    Example filter function which returns True only for prokaryotes
    """
    good_taxa = {b'Archaea', b'Bacteria', }
    taxa = re.search(b'OC.*\n', record).group()[5:]
    base_taxa = taxa.split(b'; ')[0]
    good = base_taxa in good_taxa
    return good


def parse_raw_swiss(handle, filter_fn=None, check_date=False):
    """
    Given a raw SwissProt format file containing many sequences, return an iterator of
    raw sequence strings.

    Option filter_fn argument is for a function which takes in a raw
    SwissProt format entry and returns a boolean.  If True, the string is returned
    in the iterator, if False it is not.
    """
    if not filter_fn:
        filter_fn = lambda r: True
    stream = gzip.open(handle)
    for res in _get_record(stream, check_date):
        if filter_fn(res):
            yield res
