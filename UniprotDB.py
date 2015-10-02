import collections
from cStringIO import StringIO
import blosc

import sqlalchemy
import sqlalchemy.orm
from Bio import SeqIO
from Database import UniprotProtein as UP
from Database import initialize_database


class SeqDB(collections.Mapping):

    indexes = ['id', 'genome', 'EMBL', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name', 'taxid']
    ids = ['id', 'RefSeq', 'STRING', 'GeneID', 'PIR', 'Uni_name']

    def __init__(self,
                 database='postgresql://localhost/uniprot'):
        engine = sqlalchemy.create_engine(database)
        self.Session = sqlalchemy.orm.sessionmaker(bind=engine)

    def __getitem__(self, item):
        session = self.Session()
        t = session.query(UP.raw_record).filter(
            (UP.id == item) | (UP.STRING == item) | (UP.RefSeq == item) | (UP.Uni_name == item)
        ).first()
        session.close()
        if t is None:
            return None
        r = SeqIO.read(StringIO(blosc.decompress(t[0])), 'swiss')

        return r

    def __iter__(self):
        session = self.Session()
        for row in session.query(UP.raw_record):
            yield SeqIO.read(StringIO(blosc.decompress(row[0])), 'swiss')
        session.close()

    def keys(self):
        session = self.Session()
        keys = [i[0] for i in session.query(UP.id)]
        session.close()
        return keys

    def iterkeys(self):
        session = self.Session()
        keys = (i[0] for i in session.query(UP.id))
        return keys

    def __len__(self):
        session = self.Session()
        size = session.query(UP.id).count()
        session.close()
        return size

    def get_by(self, attr, value):
        session = self.Session()
        res = session.query(UP.raw_record).filter_by(**{attr: value})
        session.close()
        return (SeqIO.read(StringIO(blosc.decompress(i[0])), 'swiss') for i in res)


def create_index(flatfiles, database='postgresql://localhost/uniprot'):
    """
    Given a list of SwissProt flatefile filenames and a SQLAlchemy database identifier,
    fill the database with the protein entires and returns a SeqDB object.
    """
    initialize_database(flatfiles, database)
    return SeqDB(database=database)

if __name__ == '__main__':
    s = SeqDB()