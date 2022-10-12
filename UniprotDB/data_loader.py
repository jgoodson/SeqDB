import itertools
import logging
import os
import tempfile
import traceback
from io import BufferedReader
from typing import Iterable, Union, BinaryIO, Tuple, Generator, List, Optional, Callable

from UniprotDB.UniprotDB import SeqDB


def open_dat(dat: str) -> Union[BufferedReader, BinaryIO]:
    """
    Opens the input string as a binary file handle depending on compression
    :param dat: String filename or url
    :return: Binary handle
    """
    if dat.startswith('http://') or dat.startswith('https://') or dat.startswith('ftp://'):
        import urllib.request
        h = urllib.request.urlopen(dat)
    else:
        h = dat
    if dat.endswith('.dat'):
        return open(h, 'rb')
    elif dat.endswith('.dat.gz') or dat.endswith('.dat.bgz'):
        import gzip
        return gzip.open(h, 'rb')
    elif dat.endswith('.dat.zst'):
        import zstandard
        import io
        cctx = zstandard.ZstdDecompressor()
        return io.BufferedReader(cctx.stream_reader(open(h, 'rb')))


def get_chunk(handle: BinaryIO, chunksize: int) -> Tuple[bytes, bool]:
    """
    Gets a chunk of data from an SwissProt format binary text handle.
    Ensures that complete records are read.
    :param handle: Binary file handle with swissprot data
    :param chunksize: Number of bytes to read before completion of record
    :return: Tuple containing (swissprot record bytes, bool whether this is last record in handle)
    """
    data = [handle.read(chunksize), handle.readline()]
    if not data[-1]:
        return data[0], True
    is_last_chunk = len(data[0]) < chunksize
    while data[-1][:3] != b'//\n':
        data.append(handle.readline())
    data = b''.join(data)
    return data, is_last_chunk


def feed_files(input_handle: BinaryIO,
               output_handles: Iterable[BinaryIO],
               chunksize: int = 50000) -> Generator[bool, None, None]:
    """
    Generator which writes reads chunks from the input handle and distributed them evenly
    among the output handles.
    :param input_handle: Binary swissprot file read handle
    :param output_handles: Iterable of binary write handles
    :param chunksize: Number of bytes to read before completion of record
    :return: bool indicating whether this is finished writing the input handle
    """
    output_cycle = itertools.cycle(output_handles)
    is_last_chunk = False
    while not is_last_chunk:
        data, is_last_chunk = get_chunk(input_handle, chunksize)
        ch = next(output_cycle)
        ch.write(data)
        yield is_last_chunk


def make_fifos(jobs: int, directory: str) -> List[str]:
    """
    Creates temporary unix FIFO named pipes
    :param jobs: Number of pipes to prepare
    :param directory: Location to make pipes
    :return: List of Filenames pointing to named pipes
    """
    fifos = []
    for i in range(jobs):
        name = os.path.join(directory, f'fifo{i}')
        os.mkfifo(name)
        fifos.append(name)
    return fifos


def process(host: str, dbtype: str, filename: str, filter_fn: Optional[Callable]=None, 
            db_kwargs=None) -> None:
    """
    Function to create a SeqDB, open a file and write the SwissProt data to the SeqDB.
    Intended for use in a multiprocessing pool
    :param host: hostname or folder location for the SeqDB
    :param dbtype: type of datastore ('lmdb', 'mongo', ...)
    :param filename: file to open for reading
    :param db_kwargs: dictionary of additional arguments for SeqDB
    :return: None
    """
    if db_kwargs is None:
        db_kwargs = {}
    s = SeqDB(host=host, dbtype=dbtype, **db_kwargs)
    try:
        with open(filename, 'rb') as fh:
            logging.debug(f"Opened {filename} for reading")
            s.update([fh], loud=False, filter_fn=filter_fn)
    except Exception as e:
        print(''.join(traceback.format_tb(e.__traceback__)))
        print(e)
        raise


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Process Uniprot flatfiles for SeqDB storage')
    parser.add_argument('dats', metavar='dats', type=str, nargs='+',
                        help='Uniprot flatfile from RefSeq (can be gzip- or zstd-compressed')
    parser.add_argument('-l', '--location', default='~/.seqdb', help='Location of the database (hostname or filename)')
    parser.add_argument('-d', '--debug', action='store_true', help='Log at DEBUG level')
    parser.add_argument('-t', '--type', default='lmdb', help='Database type to utilize')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show progress bar during storage')
    parser.add_argument('-j', '--jobs', default=8, type=int, help='Number of worker processes')
    parser.add_argument('-i', '--initialize', action='store_true', help='Whether to initialize the database')
    parser.add_argument('-n', '--num-seqs', default=0, type=int, help='Cosmetic: number of sequences in files')
    parser.add_argument('--no-index', action='store_false', help='Skip metadata indexing')
    parser.add_argument('--lmdb-db-splits', default=10, type=int, help='How many databases to split main database to')
    parser.add_argument('--lmdb-index-splits', default=10, type=int, help='How many databases to split index databases')

    args = parser.parse_args()

    logging.basicConfig(filename='data_loader.log', level=logging.DEBUG if args.debug else logging.INFO)

    process_main(dats=args.dats, 
                 location=args.location, 
                 dbtype=args.type, 
                 initialize=args.initialize, 
                 verbose=args.verbose, 
                 n_jobs=args.jobs, 
                 num_seqs=args.num_seqs,
                 db_kwargs={
                    'db_splits': args.lmdb_db_splits, 
                    'index_db_splits': args.lmdb_index_splits, 
                    'index': args.no_index
                    }
                )


def process_main(dats: Iterable[str],
                 location: str,
                 dbtype: str = 'lmdb',
                 initialize: bool = False,
                 verbose: bool = True,
                 n_jobs: int = 8,
                 num_seqs: int = 0,
                 filter_fn: Callable = lambda _: True,
                 db_kwargs: Optional[dict] = None) -> SeqDB:
    """
    Main function for parallel data loading into a SeqDB.
    :param dats: Iterable of filenames/urls for input
    :param location: hostname or folder location for the SeqDB
    :param dbtype: type of datastore ('lmdb', 'mongo', ...)
    :param initialize: bool Whether to initialize the database (remove existing data)
    :param verbose: bool Whether to show a progress bar
    :param n_jobs: number of parallel processes to use
    :param num_seqs: cosmetic number of input sequences for progress bar
    :param db_kwargs: dictionary with extra parameters for SeqDB
    :return: SeqDB object with the resulting data
    """
    from multiprocessing import get_context
    from time import sleep

    seqdb = SeqDB(host=location, dbtype=dbtype, **db_kwargs)
    if initialize:
        logging.debug('Initializing db')
        seqdb.initialize([], loud=False)

    if verbose:
        from tqdm import tqdm
        if num_seqs:
            pbar = tqdm(total=num_seqs)
        else:
            pbar = tqdm()
        current = len(seqdb)

    with tempfile.TemporaryDirectory() as directory:
        fifos = make_fifos(n_jobs, directory)

        mp_context = get_context('spawn')
        with mp_context.Pool(n_jobs) as p:
            for dat in dats:
                with open_dat(dat) as fh:
                    logging.debug('Created pool')
                    res = p.starmap_async(process, [(location, dbtype, fifo, filter_fn, db_kwargs) for fifo in fifos], chunksize=1)
                    logging.debug('Started async processing')
                    logging.debug(f'Opening fifos {fifos}')
                    output_handles = [open(f, 'wb') for f in fifos]
                    logging.debug(f'Opening fifos')
                    feed_steps = feed_files(fh, output_handles)
                    all_fed = False
                    logging.debug(f'Starting fifo feeding from {dat}')

                    while not res.ready():
                        if all_fed:
                            sleep(1)
                        else:
                            all_fed = next(feed_steps)
                            if all_fed:
                                logging.debug('Finished last write, closing FIFOs')
                                for f in output_handles:
                                    f.close()
                        if verbose:
                            new_count = len(seqdb)
                            pbar.update(new_count - current)
                            current = new_count

    return seqdb


if __name__ == '__main__':
    main()
