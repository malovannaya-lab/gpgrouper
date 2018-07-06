"""Run this directly to run multi-core grouper
Requires proper configuration of auto-grouper
Works for windows operating systems"""

import argparse
from itertools import repeat
import multiprocessing as mp
from configparser import ConfigParser
from getpass import getuser
from .pygrouper import *
from .auto_grouper import file_checker, update_database
from .cli import Config
from .containers import UserData
from . import _version


CONFIG_NAME = 'pygrouper_config.ini'
BASEDIR, _ = os.path.split(os.path.abspath(__file__))

parser = ConfigParser(comment_prefixes=(';')) # allow number sign to be read in configfile
parser.optionxform = str

def find_configfile(path='.'):
    target = os.path.join(path, CONFIG_NAME)
    if os.path.isfile(target):
        return os.path.join(target)
    else:
        return None

def get_configfile(config, config_file):
    """Get the config file for a given user"""
    if config_file is None:
        config_file = find_configfile(path=os.getcwd())
    if config_file is None:  # fail to find configfile
        return None
    config.CONFIG_FILE = config_file
    parser.read(config_file)
    return parser

def parse_configfile(config_file=None):
    """Parse the configfile and update the variables in a Config object if
    config_file exists"""
    config = Config(getuser())
    parser = get_configfile(config, config_file)
    if parser is None:
        return None
    config.inputdir = parser.get('directories', 'inputdir')
    config.outputdir = parser.get('directories', 'outputdir')
    config.rawfiledir = parser.get('directories', 'rawfiledir')
    config.contaminants = parser.get('directories', 'contaminants')

    fv_section = parser['filter values']
    filtervalues = {'ion_score': fv_section.getfloat('ion score'),
                    'qvalue': fv_section.getfloat('q value'),
                    'pep': fv_section.getfloat('PEP'),
                    'idg': fv_section.getfloat('IDG'),
                    'zmin': fv_section.getint('charge_min'),
                    'zmax': fv_section.getint('charge_max'),
                    'modi': fv_section.getint('max modis'),
                    }
    config.filtervalues = filtervalues

    column_aliases = dict()
    for column in parser['column names']:
        column_aliases[column] = [x.strip() for x in
                                  parser.get('column names', column).splitlines() if x]
    config.column_aliases = column_aliases

    refseqs = dict()
    for taxon, location in parser.items('refseq locations'):
        refseqs[int(taxon)] = location
    config.refseqs = refseqs

    labels = dict()
    for label in parser['labels']:
        labels[label] = [x.strip() for x in
                         parser.get('labels', label).splitlines() if x]
    config.labels = labels
    return config

def get_test_data():
    files = ['30259_1_rand_sample.txt', '30404_1_QEP_ML262_75min_020_RPall.txt',
             '30490_1_1_EQP_6KiP_all.txt', '30595_1_2_3T3_LM2_5_5_PROF_75MIN_all.txt']
    fdir = os.path.join(BASEDIR, 'manual_tests')
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

    current = mp.current_process()
    usrdata, databases, labels, test = params
    print('{!r} processing {!r}'.format(current, usrdata))
    sys.stdout = sys.stderr = null
    grouper(usrdata,
            database=databases[usrdata.taxonid],
            gid_ignore_file=gid_ignore_file, labels=labels)
    if test is False:
        update_database(usrdata)
    sys.stderr = _stderr
    sys.stdout = _stdout
    print('{!r} finished'.format(current, usrdata))
    return 0

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Multiprocessing pygrouper execution')
    argparser.add_argument('-p', '--processes', action="store", type=int, default=None)
    argparser.add_argument('-m', '--max-files', action="store", type=int, default=4)
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
        refseqs = {9606: os.path.join(BASEDIR, 'manual_tests/human_refseq.tab')}
    else:
        usrdatas = file_checker(INPUT_DIR, OUTPUT_DIR, maxfiles, )
    if usrdatas is None:
        print('No files')
        sys.exit(0)

    for usrdata in usrdatas:
        print('Found {!r}'.format(usrdata))

    print()
    print('\nrelease date: {}'.format(_version.__copyright__))
    print('Pygrouper v{}'.format(_version.__version__))
    print('Python version ' + sys.version)
    print('Pandas version: ' + pd.__version__)
    print('{} files found.'.format(len(usrdatas)))
    print('Running on {} processes'.format(processes))

    for usrdata in usrdatas:
        usrdata.read_csv(sep='\t')  # read from the stored psms file
        usrdata.df['metadatainfo'] = ''
        standard_names = column_identifier(usrdata.df, column_aliases)
        usrdata.df.rename(columns={v: k
                                   for k,v in standard_names.items()},
                          inplace=True)
        usrdata.df = set_modifications(usrdata.df)

    usrdatas, databases = match(usrdatas, refseqs)

    pool = mp.Pool(processes=processes)
    params = zip(usrdatas, repeat(databases), repeat(labels), repeat(gid_ignore_file), repeat(test))

    results = pool.map(work, params)


    # sys.stderr = _stderr
    # sys.stdout = _stdout
