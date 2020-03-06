import sys
from abc import ABC, abstractmethod
from functools import partial


class BaseDatabase(ABC):
    ids = ['_id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name']
    indices = ['RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'PDB', 'EMBL', 'GO', 'Pfam', 'Proteomes', 'genome',
               'taxid']

    @abstractmethod
    def __init__(self, database, host, compressor=None, decompressor=None, create_protein_func=None):
        self.database = database
        self.host = host
        if not compressor:
            import zstd
            self.compressor = zstd.ZstdCompressor()
        if not decompressor:
            import zstd
            self.decompressor = zstd.ZstdDecompressor()
        if not create_protein_func:
            from UniprotDB._utils import _create_protein_swiss
            self.create_protein_func = partial(_create_protein_swiss, compressor=self.compressor)
        from UniprotDB._utils import _extract_seqrecord
        self._extract_seqrecord = partial(_extract_seqrecord, decompressor=self.decompressor)
        pass

    def initialize(self, seq_handles, filter_fn=None, loud=False, n_seqs=None, workers=1):
        if loud:
            print("--initializating database\n", file=sys.stderr)
        self._reset()

        self._create_indices()

        self.update(seq_handles, filter_fn=filter_fn, loud=loud, total=n_seqs)

        if loud:
            print("--initialized database\n", file=sys.stderr)

    @abstractmethod
    def get_item(self, item):
        pass

    @abstractmethod
    def get_iter(self):
        pass

    @abstractmethod
    def get_iterkeys(self):
        pass

    @abstractmethod
    def get_keys(self):
        pass

    @abstractmethod
    def length(self):
        pass

    @abstractmethod
    def get_by(self, attr, value):
        pass

    @abstractmethod
    def _reset(self):
        pass

    @abstractmethod
    def _create_indices(self):
        pass

    @abstractmethod
    def update(self, handles, filter_fn=None, loud=False, total=None, processes=1):
        pass

    def add_record(self, raw_record, test=None, test_attr=None):
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
    def add_protein(self, protein):
        pass
