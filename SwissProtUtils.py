from __future__ import print_function, division
from sys import stderr
from Bio import bgzf


def _get_record(handle):
    """
    Returns the next complete SwissProt entry in the input handle
    """
    lines = []
    for line in handle:
        lines.append(line)
        if line.startswith('//'):
            return ''.join(lines)
        if line == '':
            return None


def filter_proks(record):
    """
    Example filter function which returns True only for prokaryotes
    """
    good_taxa = {'Archaea', 'Bacteria', }
    for l in record.split('\n'):
        if l.startswith('OC'):
            taxa = l[5:]
            break
    base_taxa = taxa.split('; ')[0]
    good = base_taxa in good_taxa
    return good


def parse_raw_swiss(filename, filter_fn=None):
    """
    Given a raw SwissProt format file containing many sequences, return an iterator of
    raw sequence strings.

    Option filter_fn argument is for a function which takes in a raw
    SwissProt format entry and returns a boolean.  If True, the string is returned
    in the iterator, if False it is not.
    """
    if not filter_fn:
        filter_fn = lambda r: True
    handle = bgzf.open(filename)
    while True:
        res = _get_record(handle)
        if not res:
            break
        if filter_fn(res):
            yield res

if __name__ == '__main__':
    it = parse_raw_swiss('/Volumes/OldSynoRaid/protein/uniprot_sprot.prok.dat.bgz', filter_fn=filter_proks)
    print(it.next())