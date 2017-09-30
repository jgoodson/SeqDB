from __future__ import print_function, division

import gzip
import re


def _get_record(handle):
    """
    Returns the next complete SwissProt entry in the input handle

    :param handle: file handle
    :return: string containing single SwissProt entry
    """
    lines = []
    for line in handle:
        lines.append(line)
        if line.startswith(b'//'):
            return b''.join(lines)
        if line == b'':
            return None


def filter_proks(record):
    """
    Example filter function which returns True only for prokaryotes

    :param record: string containing single SwissProt entry
    :return: bool result
    """
    good_taxa = {b'Archaea', b'Bacteria', }
    taxa = re.search(b'OC.*\n', record).group()[5:]
    base_taxa = taxa.split(b'; ')[0]
    good = base_taxa in good_taxa
    return good


def parse_raw_swiss(filename, filter_fn=None):
    """
    Given a raw SwissProt format file containing many sequences, return an iterator of
    raw sequence strings.

    Option filter_fn argument is for a function which takes in a raw
    SwissProt format entry and returns a boolean.  If True, the string is returned
    in the iterator, if False it is not.

    :param filename: string filename for Uniprot flatfile
    :param filter_fn: function that returns True/False given flatfile record
    :return: iterator with individual SwissProt records
    """
    if not filter_fn:
        def filter_fn(*args):
            return True
    # handle = gzip.open(filename)

    with gzip.open(filename) as handle:
        while True:
            res = _get_record(handle)
            if not res:
                break
            if filter_fn(res):
                yield res
