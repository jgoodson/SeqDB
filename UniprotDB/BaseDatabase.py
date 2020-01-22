from abc import ABC, abstractmethod
import sys

class BaseDatabase(ABC):
    ids = ['_id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name']

    @abstractmethod
    def __init__(self, database, host):
        pass

    def initialize(self, seq_handles, filter_fn=None, loud=False, n_seqs=None, processes=1):
        if loud:
            print("--initializating database\n", file=sys.stderr)
        self._reset()

        self.update(seq_handles, filter_fn=filter_fn, total=n_seqs, processes=processes)

        self._create_indices()

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

    @abstractmethod
    def add_protein(self, protein, test=None, test_attr=None):
        pass