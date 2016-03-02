import os
from configparser import ConfigParser
import click
from pygrouper.manual_tests import test as manual_test
from pygrouper import subfuncts
from pygrouper import auto_grouper
__author__ = 'Alexander B. Saltzman'
__copyright__ = 'Copyright January 2016'
__credits__ = ['Alexander B. Saltzman', 'Anna Malovannaya']
__license__ = 'MIT'
__version__ = '0.1.015'
__maintainer__ = 'Alexander B. Saltzman'
__email__ = 'saltzman@bcm.edu'


#HOMEDIR = os.path.expanduser('~')
PROFILE_DIR = click.get_app_dir('pygrouper', roaming=False, force_posix=True)
#PROFILE_DIR = os.path.join(HOMEDIR, '.pygrouper')
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

class Config(object):

    def __init__(self, user='profile_default'):
        self.user = user
        self.ispec_url = None
        self.database = None
        self.outfile = '-'
        self.inputdir = PROFILE_DIR
        self.filtervalues = dict()
        self.column_aliases = dict()
        self.CONFIG_DIR = PROFILE_DIR
        self.inputdir = PROFILE_DIR
        self.outputdir = PROFILE_DIR
        self.rawfiledir = PROFILE_DIR
        self.refseqs = dict()
class CaseConfigParser(ConfigParser):
    def optionxform(self, optionstr):
        return optionstr
parser = ConfigParser(comment_prefixes=(';')) # allow number sign to be read in configfile
parser.optionxform = str


@click.group()
@click.version_option(__version__)
@click.option('-p', '--profile', type=str, default='profile_default',
              help='Select the profile to use')
@click.pass_context
def cli(ctx, profile):

    if ctx.invoked_subcommand is not None:
        ctx.obj = Config(profile)
        parse_configfile()

pass_config = click.make_pass_decorator(Config)

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument('profile', type=str)
@pass_config
def makeuser(config, profile):
    """Make a new profile with new config file.
    Note this happens automatically when chaining the
    --profile option with any subcommand."""
    config.user = profile
    #click.echo('making user')
    parse_configfile()


@pass_config
def make_configfile(config):
    """Generate a new config file
    (and also parent directory if necessary)"""
    if not os.path.isdir(PROFILE_DIR):
        os.mkdir(PROFILE_DIR)
    CONFIG_DIR = os.path.join(PROFILE_DIR, config.user)
    config.CONFIG_DIR = CONFIG_DIR
    if not os.path.isdir(CONFIG_DIR):
        os.mkdir(CONFIG_DIR)
    BASE_CONFIG = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               'base_config.ini')
    #BASE_CONFIG = os.path.join()
    click.echo(BASE_CONFIG)
    parser.read(BASE_CONFIG)
    parser.set('profile', 'user', config.user)
    CONFIG_FILE = os.path.join(CONFIG_DIR, 'config.ini')
    #click.echo(CONFIG_FILE)
    write_configfile(CONFIG_FILE, parser)
    click.echo('Creating new profile : {}'.format(config.user))

def write_configfile(CONFIG_FILE, parser):
    """Write/update config file with parser"""
    #CONFIG_DIR = config.CONFIG_DIR
    with open(CONFIG_FILE, 'w') as f:
        parser.write(f)

@pass_config
def get_configfile(config):
    """Get the config file for a given user"""
    CONFIG_DIR = os.path.join(PROFILE_DIR, config.user)
    CONFIG_FILE = os.path.join(CONFIG_DIR, 'config.ini')
    config.CONFIG_DIR = CONFIG_DIR
    config.CONFIG_FILE = CONFIG_FILE
    if not os.path.isfile(CONFIG_FILE):
        make_configfile()
    parser.read(CONFIG_FILE)
    return parser

@pass_config
def parse_configfile(config):
    """Parse the configfile and update the variables in a Config object"""
    parser = get_configfile()
    config.inputdir = parser.get('directories', 'inputdir')
    config.outputdir = parser.get('directories', 'outputdir')
    config.rawfiledir = parser.get('directories', 'rawfiledir')
    fv_section = parser['filter values']
    filtervalues = {'Filter_IS': fv_section.getfloat('ion score'),
                    'Filter_qV': fv_section.getfloat('q value'),
                    'Filter_PEP': fv_section.getfloat('PEP'),
                    'Filter_IDG': fv_section.getfloat('IDG'),
                    'Filter_Z_min': fv_section.getint('charge_min'),
                    'Filter_Z_max': fv_section.getint('charge_max'),
                    'Filter_Modi': fv_section.getint('max modis'),
                    }
    config.filtervalues = filtervalues
    column_aliases = dict()
    for column in parser['column names']:
        column_aliases[column] = [x.strip() for x in
                               parser.get('column names', column).splitlines() if x]
    config.column_aliases = column_aliases
    refseqs = dict()
    for taxon, location in parser.items('refseq locations'):
        refseqs[taxon] = {'loc': location,
                          'size': parser.getfloat('refseq file sizes', taxon)} 
    config.refseqs = refseqs


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-q', '--quick', is_flag=True, help='Run a single small file')
@click.option('-p', '--profile', is_flag=True, help='Run a large profiling file with multiple taxons')
@pass_config
def test(config, quick, profile ):
    """Test pygrouper with some pre-existing data."""
    parse_configfile()
    INPUT_DIR = config.inputdir
    OUTPUT_DIR = config.outputdir
    RAWFILE_DIR = config.rawfiledir
    refseqs = config.refseqs
    filtervalues = config.filtervalues
    column_aliases = config.column_aliases
    gid_ignore_file = os.path.join(config.CONFIG_DIR, 'geneignore.txt')
    manual_test.runtest(quick, profile, inputdir=INPUT_DIR, outputdir=OUTPUT_DIR,
                        rawfilepath=RAWFILE_DIR, refs=refseqs, FilterValues=filtervalues,
                        column_aliases=column_aliases, gid_ignore_file=gid_ignore_file,
                        configpassed=True)

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument('taxonid', type=int)
@click.argument('taxonfile', type=click.Path(exists=True))
@pass_config
def add_taxon(config, taxonid, taxonfile):
    """Add a taxon to for use with pygrouper.
    Sets the taxonid -> taxonfile which is the file to use for the given taxonid."""

    filesize = str(subfuncts.bufcount(taxonfile))
    parser = get_configfile()
    parser.set('refseq locations', str(taxonid), taxonfile)
    parser.set('refseq file sizes', str(taxonid), filesize)
    write_configfile(config.CONFIG_FILE, parser)
    click.echo('Updating taxon id {} to path {}'.format(taxonid, taxonfile))

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('--path-type', type=click.Choice(['input', 'output', 'rawfile']),
              default='outputdir')
@click.argument('path', type=click.Path(exists=True))
@pass_config
def setpath(config, path_type, path):
    """Set a path for pygrouper"""
    parser = get_configfile()
    category = path_type+'dir'
    parser.set('directories', category, path)
    write_configfile(config.CONFIG_FILE, parser)
    click.echo('Updating {} to {}'.format(category, path))

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-a', '--autorun', is_flag=True,
              help='Run automatically by scanning directory and connecting to iSPEC')
@click.option('-m', '--max-files', type=int, default=99,
              help='Maximum number of experiments to quene for autorun')
@pass_config
def run(config, autorun, max_files):
    """Run PyGrouper"""
    parse_configfile()
    INPUT_DIR = config.inputdir
    OUTPUT_DIR = config.outputdir
    RAWFILE_DIR = config.rawfiledir
    refseqs = config.refseqs
    filtervalues = config.filtervalues
    column_aliases = config.column_aliases
    gid_ignore_file = os.path.join(config.CONFIG_DIR, 'geneignore.txt')
    if autorun:
        auto_grouper.interval_check(3600, INPUT_DIR, OUTPUT_DIR,
                                    max_files, rawfilepath=RAWFILE_DIR,
                                    refs=refseqs, FilterValues=filtervalues,
                                    column_aliases=column_aliases,
                                    gid_ignore_file=gid_ignore_file,
                                    configpassed=True)
    else:
        click.echo('Not implemented yet')
        usrfile = click.prompt('Enter a file to group', type=click.Path(exists=True))
