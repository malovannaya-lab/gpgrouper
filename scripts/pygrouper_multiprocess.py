"""Run this directly to run multi-core grouper
Requires proper configuration of auto-grouper
Works for windows operating systems"""

import os
import sys
import argparse
from datetime import datetime
from itertools import repeat
from configparser import ConfigParser
import multiprocessing as mp

import pandas as pd

from pygrouper import pygrouper
from pygrouper.auto_grouper import file_checker, update_database
from pygrouper.cli import Config
from pygrouper.containers import UserData
from pygrouper.parse_config import parse_configfile
from pygrouper import _version


# BASEDIR, _ = os.path.join(os.path.dirname(__file__), '..', 'manual_tests')
BASEDIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', 'pygrouper/manual_tests')
BASEDIR = os.path.abspath(BASEDIR)

parser = ConfigParser(comment_prefixes=(';')) # allow number sign to be read in configfile
parser.optionxform = str


def get_test_data():
    files = ['30259_1_rand_sample.txt',
             # '30404_1_QEP_ML262_75min_020_RPall.txt',
             # '30490_1_1_EQP_6KiP_all.txt', '30595_1_2_3T3_LM2_5_5_PROF_75MIN_all.txt'
    ]
    # fdir = os.path.join(BASEDIR, 'manual_tests')
    fdir = BASEDIR
    if not os.path.exists(os.path.join(fdir, 'out')):
        os.mkdir(os.path.join(fdir, 'out'))

    usrdatas = list()
    for f in files:
        usrdata = UserData()
        usrdata.recno = f[0:5]
        usrdata.outdir = os.path.join(fdir, 'out')
        usrdata.quant_source = 'AUC'
        usrdata.filtervalues['ion_score'] = 7
        usrdata.filtervalues['qvalue']    = 0.05
        usrdata.filtervalues['pep']       = 'all'
        usrdata.filtervalues['idg']       = 'all'
        usrdata.filtervalues['zmin']      = 2
        usrdata.filtervalues['zmax']      = 6
        usrdata.filtervalues['modi']      = 4
        usrdata.searchdb = os.path.join(fdir, 'human_refseq.tab')
        usrdata.datafile = os.path.join(fdir, f)
        usrdata.taxonid = 9606
        usrdatas.append(usrdata)
    return usrdatas


def work(params):

    null = open(os.devnull,'w')
    _stderr = sys.stderr
    _stdout = sys.stdout

    current = mp.current_process().name
    usrdata, databases, labels, gid_ignore_file, test = params
    print('{} | {} processing {!r}'.format(datetime.now(), current, usrdata))
    sys.stdout = sys.stderr = null
    error = 0
    try:
        pygrouper.grouper(usrdata,
                          database=databases[usrdata.taxonid],
                          gid_ignore_file=gid_ignore_file, labels=labels)
    except Exception as e:
        error = e

    sys.stderr = _stderr
    sys.stdout = _stdout
    print('{} finished'.format(current, usrdata))
    return {repr(usrdata): error}

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Multiprocessing pygrouper execution')
    argparser.add_argument('-p', '--processes', action="store", type=int, default=None)
    argparser.add_argument('-m', '--max-files', action="store", type=int, default=999)
    argparser.add_argument('-t', '--test', action="store_true")
    args = argparser.parse_args()
    maxfiles = args.max_files
    processes = args.processes or mp.cpu_count() - 1
    test = args.test

    config = parse_configfile(None)

    INPUT_DIR = config.inputdir or '.'
    OUTPUT_DIR = config.outputdir or '.'
    RAWFILE_DIR = config.rawfiledir or '.'
    LABELS = config.labels
    refseqs = config.refseqs
    filtervalues = config.filtervalues
    column_aliases = config.column_aliases
    gid_ignore_file = config.contaminants

    if test:
        usrdatas = get_test_data()
        refseqs = {9606: os.path.join(BASEDIR, 'human_refseq.tab')}
    else:
        usrdatas = file_checker(INPUT_DIR, OUTPUT_DIR, maxfiles, )
    if usrdatas is None:
        print('No files')
        sys.exit(0)

    for usrdata in usrdatas:
        print('Found {!r}'.format(usrdata))

    startTime = datetime.now()
    print()
    print('Pygrouper v{}'.format(_version.__version__))
    print('\nrelease date: {}'.format(_version.__copyright__))
    print('Python version ' + sys.version)
    print('Pandas version: ' + pd.__version__)
    print('{} files found.'.format(len(usrdatas)))
    print('Running on {} processes'.format(processes))
    print('\nStart at {}'.format(startTime))

    usrdatas = pygrouper.set_up(usrdatas, column_aliases)

    usrdatas, databases = pygrouper.match(usrdatas, refseqs)

    pool = mp.Pool(processes=processes)
    params = zip(usrdatas, repeat(databases), repeat(LABELS), repeat(gid_ignore_file), repeat(test))


    results = pool.map(work, params)
    errors = dict()
    if test:
        sys.exit(0)

    for d in results:
        errors.update(d)
    for usrdata in usrdatas:
        if errors[repr(usrdata)] == 0:
            update_database(usrdata)
        else:
            print('Failure for {!r}'.format(usrdata))
            print(errors[repr(usrdata)])
            print()
    print('Time taken : {}\n'.format(datetime.now() - startTime))
