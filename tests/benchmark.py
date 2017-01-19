import sys
import os
import time
from io import StringIO
from collections import Counter
import traceback
import unittest

import numpy as np
import pandas as pd
from pyteomics import parser
import click
from click.testing import CliRunner
from pygrouper import cli
np.random.seed(0)

BASEDIR, _ = os.path.split(os.path.abspath(__file__))
CONFIG_FILE = os.path.join(BASEDIR, '../pygrouper_config.ini')

# from http://www.tiem.utk.edu/~gross/bioed/webmodules/aminoacid.htm
AA_FREQ = dict(
    A = 0.074,
    R = 0.042,
    N = 0.044,
    D = 0.059,
    C = 0.033,
    E = 0.058,
    Q = 0.037,
    G = 0.074,
    H = 0.029,
    I = 0.038,
    L = 0.076,
    K = 0.072,
    M = 0.018,
    F = 0.040,
    P = 0.050,
    S = 0.081,
    T = 0.062,
    W = 0.013,
    Y = 0.032,  # rounded down so sum is exactly equal to 1
    V = 0.068,
)

AAs = list(AA_FREQ.keys())
AA_prob = list(AA_FREQ.values())

def timed(n=30):
    '''
    Running a microbenchmark. Never use this.
    '''
    def deco(func):
        def wrapper(*args, **kwargs):
            timings = []
            for i in range(n):
                t0 = time.time()
                func(*args, **kwargs)
                t1 = time.time()
                timings.append(t1 - t0)
            print('for {} runs, mean : {}, std : {}'.format(n,
                                                            np.mean(timings),
                                                            np.std(timings)))
            return timings
        return wrapper
    return deco


def get_peptide(length):
    return ''.join(np.random.choice(AAs, size=int(length), p=AA_prob))

def trypsin(aa):
    return parser.cleave(aa, parser.expasy_rules['trypsin'], min_length=7, missed_cleavages=1)

def make_database(n=200, seed=None):
    database_size = n
    gids    = np.random.randint(2, 20000, size=database_size)
    protgis = np.random.randint(200000, 800000, size=database_size)

    lb = 5
    lengths_ = np.random.randn(database_size*5)
    pos_lengths = np.multiply(np.where(lengths_ > 0,
                                       lengths_, np.abs(lengths_)), 30)
    adjust = lb - np.min(pos_lengths)
    # approximately normal distribution, though in reality it's right tailed
    lengths = (pos_lengths + adjust).round()

    pept_pool = [  get_peptide(length) for length in lengths  ]
    seqs = [ ''.join(np.random.choice(pept_pool,
                                      np.random.randint(5,12)))
                     for _ in range(database_size)]

    data = dict(TaxonID = 9606,
                HomologeneID = '',
                GeneID    = gids,
                ProteinGI = protgis,
                FASTA     = seqs
    )
    return pd.DataFrame(data).sort_values(by='FASTA')

def make_psms(ref_db, n, seed=None):
    seqs_ = set()
    for s in ref_db['FASTA']:
        seqs_ |= trypsin(s)
    valid_seqs = [x for x in seqs_ if len(x) < 50]
    seqs = np.random.choice(valid_seqs, size=n, replace=True)
    areas = np.random.randint(1, 100, size=n)
    ionscores_ = np.random.normal(loc=30, scale=10, size=n)
    # approximately normal distribution centered around 30. Importantly no zeros
    ionscores = np.where(ionscores_ > 0, ionscores_, np.abs(ionscores_))
    qvalue_ = np.random.gamma(shape=1, scale=2, size=n)
    maxq = np.max(qvalue_)
    qvalue = np.divide(qvalue_, maxq) * .05

    d = dict(Sequence       = seqs,
             Modifications  = ' ',
             PrecursorArea  = np.multiply(areas, np.power(10, 7)),
             Charge         = np.random.choice([2, 3, 4], n),
             IonScore       = ionscores,
             q_value        = qvalue,
             PEP            = qvalue * 10,  # rough approximation
             SpectrumFile   = 'File1.raw',
             RTmin          = 99 * np.random.random(n) + 1,
             DeltaMassPPM   = np.random.normal(scale=.5, size=n),
             )

    return pd.DataFrame(d)

class RunThrough():

    stdout = sys.stdout
    stderr = sys.stderr
    DB  = './testdata/random_db.tab'
    FILE = './testdata/10101_1_1_random.tab'

    def __init__(self, N, seed=None, create=True, delete=True):
        self.N = N
        if seed is not None:
            self.seed = int(seed)
        else:
            self.seed = None
        self.create = create
        self.delete = delete
    def setUp(self):
        if not os.path.exists('./testdata'):
            os.mkdir('./testdata')
        db = make_database(seed=self.seed)
        psms = make_psms(db, self.N, seed=self.seed)
        print('Writing', self.DB)
        db.to_csv(self.DB, sep='\t', index=False)
        print('Writing', self.FILE)
        psms.to_csv(self.FILE, sep='\t', index=False)

    def tearDown(self):
        print('Cleaning up ...', end='')
        os.remove(self.DB)
        os.remove(self.FILE)
        for f in os.listdir('.'):
            if f.startswith('10101_1_1'):
                os.remove(f)
        for f in os.listdir('./testdata'):
            if f.startswith('10101_1_1'):
                os.remove(os.path.join('./testdata', f))
        print('done')
    def run(self):
        # if self.create or not (os.path.exists(self.DB) and
        #                        os.path.exists(self.FILE)):
        # self.setUp()
        self._exec()
        # if self.delete:
        #     self.tearDown()
        return 0

    @timed(3)
    def _exec(self):
        sys.stdout = StringIO()
        sys.stderr = StringIO()
        print('hi')
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', self.DB,
                                           '--psms-file', self.FILE,
                                           '--taxonid', 9606,
                                           '--outdir', './testdata',
                                           # '--configfile', CONFIG_FILE,
                                           '--labeltype', 'none',
                                           ],
                                 catch_exceptions=False)
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        try:
            assert 0 == response.exit_code
        except AssertionError:
            raise AssertionError('\n{}\n{!r}'.format(''.join(traceback.format_tb(response.exc_info[-1])),
                                                     response.exc_info[1])
            )

@click.command()
@click.option('--create/--no-create', default=False, show_default=True,
              help="""Create a new database.
              Will happen automatically if does not exist.""")
@click.option('--delete/--no-delete', default=True, show_default=True,
              help='Delete all output files.')
@click.option('-n', '--number', default=1000, show_default=True)
@click.option('-s', '--seed', default=None, show_default=True)
def main(create, delete, number, seed):
    RunThrough(number, seed, create=create,
               delete=delete).run()

if __name__ == '__main__':
    main()
