import gzip
import re
from datetime import datetime


def _get_record(handle, ignore={b'R', b'C'}):
    """
    Returns the next complete SwissProt entry in the input handle
    """
    lines = []
    for line in handle:
        if not line[0] in ignore:
            lines.append(line)
        if line.startswith(b'//'):
            yield b''.join(lines)
            lines = []


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
    for res in _get_record(stream):
        if filter_fn(res):
            yield res
