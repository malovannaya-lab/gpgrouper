import os

import six
if six.PY3:
    from configparser import ConfigParser
elif six.PY2:
    from ConfigParser import ConfigParser
from getpass import getuser


CONFIG_NAME = 'pygrouper_config.ini'

class Config(object):

    def __init__(self, user):
        self.user = user
        self.ispec_url = None
        self.database = None
        self.outfile = '-'
        self.filtervalues = dict()
        self.column_aliases = dict()
        self.CONFIG_DIR = '.'
        self.inputdir = '.'
        self.outputdir = '.'
        self.rawfiledir = '.'
        self.labels = dict()
        self.refseqs = dict()
        self.fastadb = None
        self.contaminants = None

if six.PY3:
    parser = ConfigParser(comment_prefixes=(';')) # allow number sign to be read in configfile
elif six.PY2:
    parser = ConfigParser() # allow number sign to be read in configfile
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

def parse_configfile(config_file=None, user=None):
    """Parse the configfile and update the variables in a Config object if
    config_file exists"""
    user = user or getuser()
    config = Config(getuser())
    parser = get_configfile(config, config_file)
    if parser is None:
        return config
    config.inputdir = parser.get('directories', 'inputdir')
    config.outputdir = parser.get('directories', 'outputdir')
    config.rawfiledir = parser.get('directories', 'rawfiledir')
    config.contaminants = parser.get('directories', 'contaminants')

    # fv_section = parser['filter values']
    # filtervalues = {'ion_score': fv_section.getfloat('ion score'),
    #                 'qvalue': fv_section.getfloat('q value'),
    #                 'pep': fv_section.getfloat('PEP'),
    #                 'idg': fv_section.getfloat('IDG'),
    #                 'zmin': fv_section.getint('charge_min'),
    #                 'zmax': fv_section.getint('charge_max'),
    #                 'modi': fv_section.getint('max modis'),
    #                 }
    filtervalues = {'ion_score': parser.getfloat('filter values', 'ion score'),
                    'qvalue':    parser.getfloat('filter values', 'q value'),
                    'pep':       parser.getfloat('filter values', 'PEP'),
                    'idg':       parser.getfloat('filter values', 'IDG'),
                    'zmin':      parser.getint('filter values', 'charge_min'),
                    'zmax':      parser.getint('filter values', 'charge_max'),
                    'modi':      parser.getint('filter values', 'max modis'),
                    }
    config.filtervalues = filtervalues

    column_aliases = dict()
    # for column in parser['column names']:
    for column, _ in parser.items('column names'):
        column_aliases[column] = [x.strip() for x in
                                  parser.get('column names', column).splitlines() if x]
    config.column_aliases = column_aliases

    refseqs = dict()
    for taxon, location in parser.items('refseq locations'):
        refseqs[int(taxon)] = location
    config.refseqs = refseqs

    labels = dict()
    for label, _ in parser.items('labels'):
        labels[label] = [x.strip() for x in
                         parser.get('labels', label).splitlines() if x]
    config.labels = labels
    return config
