import sys
from abc import ABC, abstractmethod
from functools import partial
from typing import Union, Callable, Iterable, Generator, List

import zstd
from Bio.SeqRecord import SeqRecord


class BaseDatabase(ABC):
    ids = ['_id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name']
    indices = ['RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'PDB', 'EMBL', 'GO', 'Pfam', 'Proteomes', 'genome',
               'taxid']

    @abstractmethod
    def __init__(self, database: str, host: Union[tuple, str],
                 compressor: zstd.ZstdCompressor = None,
                 decompressor: zstd.ZstdDecompressor = None,
                 create_protein_func: Callable = None):
        self.database = database
        self.host = host
        if not compressor:
            self.compressor = zstd.ZstdCompressor()
        if not decompressor:
            self.decompressor = zstd.ZstdDecompressor()
        if not create_protein_func:
            from UniprotDB._utils import _create_protein_swiss
            self.create_protein_func = partial(_create_protein_swiss, compressor=self.compressor)
        else:
            self.create_protein_func = partial(create_protein_func, compressor=self.compressor)
        from UniprotDB._utils import _extract_seqrecord
        self._extract_seqrecord = partial(_extract_seqrecord, decompressor=self.decompressor)
        pass

    def initialize(self, seq_handles: Iterable,
                   filter_fn: Callable[[bytes], bool] = None,
                   loud: bool = False,
                   n_seqs: int = None,
                   workers: int = 1) -> None:
        if loud:
            print("--initializating database\n", file=sys.stderr)
        self._reset()

        self._create_indices()

        self.update(seq_handles, filter_fn=filter_fn, loud=loud, total=n_seqs)

        if loud:
            print("--initialized database\n", file=sys.stderr)

    @abstractmethod
    def get_item(self, item: str) -> SeqRecord:
        pass

    @abstractmethod
    def get_iter(self) -> Generator[SeqRecord, None, None]:
        pass

    @abstractmethod
    def get_iterkeys(self) -> Generator[str, None, None]:
        pass

    @abstractmethod
    def get_keys(self) -> List[str]:
        pass

    @abstractmethod
    def length(self) -> int:
        pass

    @abstractmethod
    def get_by(self, attr: str, value: str) -> List[SeqRecord]:
        pass

    @abstractmethod
    def _reset(self) -> None:
        pass

    @abstractmethod
    def _create_indices(self) -> None:
        pass

    @abstractmethod
    def update(self, handles: Iterable,
               filter_fn: Callable[[bytes], bool] = None,
               loud: bool = False,
               total: int = None,
               workers: int = 1) -> None:
        pass

    def add_record(self, raw_record: bytes, test: str = None, test_attr: str = None) -> bool:
        protein = self.create_protein_func(raw_record)
        if test:
            good = False
            if test == protein['_id']:
                good = True
            if not good:
                for ref in ([test_attr] if test_attr else self.ids):
                    if test in protein.get(ref, []):
                        good = True
            if not good:
                return False
        self.add_protein(protein)
        return True

    @abstractmethod
    def add_protein(self, protein: dict) -> None:
        pass
