import re
from typing import BinaryIO, Collection, Callable


def _get_record(handle: BinaryIO, ignore: Collection[bytes] = (b'R', b'C')):
    """
    Returns the next complete SwissProt entry in the input handle
    """
    lines = []
    for line in handle:
        if line[0] not in ignore:
            lines.append(line)
        if line.startswith(b'//'):
            yield b''.join(lines)
            lines = []


def filter_proks(record: bytes):
    """
    Example filter function which returns True only for prokaryotes
    """
    good_taxa = {b'Archaea', b'Bacteria', }
    taxa = re.search(b'OC.*\n', record).group()[5:]
    base_taxa = taxa.split(b'; ')[0]
    return base_taxa in good_taxa


def parse_raw_swiss(handle: BinaryIO, filter_fn: Callable[[bytes], bool] = None):
    """
    Given a raw SwissProt format file containing many sequences, return an iterator of
    raw sequence strings.

    Option filter_fn argument is for a function which takes in a raw
    SwissProt format entry and returns a boolean.  If True, the string is returned
    in the iterator, if False it is not.
    """
    if not filter_fn:
        def filter_fn(*args):
            return True
    for res in _get_record(handle):
        if filter_fn(res):
            yield res
