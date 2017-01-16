import os
import sys
import re
from itertools import zip_longest
from getpass import getuser
from pathlib import Path
from datetime import datetime
from configparser import ConfigParser
import click

# from .manual_tests import test as manual_test
from . import subfuncts, auto_grouper, pygrouper, _version
from pygrouper import genericdata as gd
from pygrouper.containers import UserData
from .parse_config import parse_configfile, find_configfile, Config


# from manual_tests import test as manual_test
# import subfuncts, auto_grouper, grouper
# from pygrouper import genericdata as gd
__author__ = 'Alexander B. Saltzman'
__copyright__ = 'Copyright January 2016'
__credits__ = ['Alexander B. Saltzman', 'Anna Malovannaya']
__license__ = 'MIT'
__version__ = _version.__version__
__maintainer__ = 'Alexander B. Saltzman'
__email__ = 'saltzman@bcm.edu'


#HOMEDIR = os.path.expanduser('~')
CONFIG_DIR = click.get_app_dir('pygrouper', roaming=False, force_posix=True)
#PROFILE_DIR = os.path.join(HOMEDIR, '.pygrouper')
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
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

parser = ConfigParser(comment_prefixes=(';')) # allow number sign to be read in configfile
parser.optionxform = str


@click.group(name='main')
@click.version_option(__version__)
@click.option('-p', '--profile', type=str, default=getuser(),
              help='Name of the user.')
@click.pass_context
def cli(ctx, profile):
    pass


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-p', '--path', type=click.Path(exists=True), default='.')
def make_configfile(path):
    """Generate a new config file
    (and also parent directory if necessary)"""
    if not os.path.isdir(CONFIG_DIR):
        os.mkdir(CONFIG_DIR)
    # CONFIG_DIR = os.path.join(PROFILE_DIR, config.user)
    # config.CONFIG_DIR = CONFIG_DIR
    # if not os.path.isdir(CONFIG_DIR):
        # os.mkdir(CONFIG_DIR)
    BASE_CONFIG = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               'base_config.ini')
    parser.read(BASE_CONFIG)
    parser.set('profile', 'user', config.user)
    CONFIG_FILE = os.path.join(CONFIG_DIR, 'pygrouper_config.ini')

    write_configfile(os.path.join(path, 'pygrouper_config.ini'), parser)
    click.echo('Creating new config file.')

def write_configfile(CONFIG_FILE, parser):
    """Write/update config file with parser"""
    #CONFIG_DIR = config.CONFIG_DIR
    with open(CONFIG_FILE, 'w') as f:
        parser.write(f)

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-p', '--path', type=click.Path(exists=True), default='.')
def openconfig(path):
    config_file = find_configfile(path=path)
    if config_file is None:
        click.echo("Could not find config file in {}".format(path))
        return
    click.launch(CONFIG_FILE)

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-q', '--quick', is_flag=True, help='Run a single small file')
@click.option('-p', '--profile', is_flag=True, help='Run a large profiling file with multiple taxons')
@click.option('-t', '--tmt', is_flag=True, help='Run a large profiling file with multiple taxons')
@click.pass_obj
def test(config, quick, profile, tmt):
    """Test pygrouper with some pre-existing data."""
    parse_configfile()
    INPUT_DIR = config.inputdir
    OUTPUT_DIR = config.outputdir
    RAWFILE_DIR = config.rawfiledir
    LABELS = config.labels
    refseqs = config.refseqs
    filtervalues = config.filtervalues
    column_aliases = config.column_aliases
    gid_ignore_file = os.path.join(config.CONFIG_DIR, 'geneignore.txt')
    manual_test.runtest(quick, profile, tmt, inputdir=INPUT_DIR, outputdir=OUTPUT_DIR,
                        rawfilepath=RAWFILE_DIR, refs=refseqs, FilterValues=filtervalues,
                        column_aliases=column_aliases, gid_ignore_file=gid_ignore_file,
                        labels=LABELS,)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-p', '--path', type=click.Path(exists=True), default='.')
@click.argument('taxon', type=int, nargs=-1)
def download_taxon(path, taxon):
    """Downloads a new taxon"""
    if any(x not in (10090, 9606) for x in taxon):
        raise NotImplementedError('No support for updating taxon {} yet.'.format(*(x for x in taxon
                                                                                   if x not in (10090,
                                                                                                9606))))
    gz_path = os.path.join(CONFIG_DIR, 'tempfiles')
    print(path)
    if not os.path.exists(gz_path):
        os.mkdir(gz_path)

    for download in (gd.download_ebi_files, gd.download_ncbi_files):
        try:
            download(path=gz_path, taxa=taxon)
        except Exception as e:
            # gd.cleanup(gz_path)
            raise(e)
    gd.unzip_all(path=gz_path)
    gd.entrylist_formatter(path=gz_path)
    gd.protein2ipr_formatter(path=gz_path)
    for taxa in taxon:
        gd.idmapping_formatter(taxa, path=gz_path)
        gd.inputfiles = append_all_files(taxa, path=gz_path)
        gd.refseq = refseq_dict(inputfiles) # Make refseq dictionary
        gd.gene2accession_formatter(taxa, path=gz_path)
        g2a = gd.g2acc_dict(refseq, path=gz_path)
        gd.homologene_formatter(path=gz_path)
        hid = gd.hid_dict(path=gz_path)
        gd.file_input(inputfiles, refseq, g2a, hid, taxonid)
        gd.file_write(taxa, lines_seen, path=path)
        gd.refseq_formatter_4_mascot(taxon, path=path)

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-p', '--path', type=click.Path(exists=True), default='.')
def view_taxons(path):
    """List current taxa and their locations based on the config file"""
    parser = get_configfile(path)
    for taxon, file in parser.items('refseq locations'):
        filestat = os.stat(file)
        click.echo('{} : {}\t{}'.format(taxon, file,
                                        datetime.fromtimestamp(filestat.st_mtime)))

DEFAULTS = {'max_files': 99, 'pep': 1.0, 'enzyme': 'trypsin', 'configfile': None, 'interval': 3600, 'rawfiledir': '.',
            'taxonid': None, 'contaminants': None, 'quant_source': 'AUC', 'outdir': None, 'zmax': 6, 'ion_score': 7.0,
            'ion_score_bins': (10.0, 20.0, 30.0),
            'qvalue': 0.05, 'idg': 9, 'autorun': False, 'name': 'shiro', 'modi': 4, 'zmin': 2,
            'no_taxa_redistrib': False, 'labeltype': 'none', 'pipeline': 'PD'}

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-a', '--autorun', is_flag=True,
              help='Run automatically by scanning directory and connecting to iSPEC database.')
@click.option('-c', '--contaminants', type=click.Path(exists=True, dir_okay=False),
              help='Contaminants file of IDs to ignore when calculating with multiple taxa.')
@click.option('-d', '--database', type=click.Path(exists=True, dir_okay=False),
              help='Database file to use. Ignored with autorun.')
@click.option('-e', '--enzyme', type=click.Choice(['trypsin', 'trypsin/p', 'chymotrypsin',
                                                   'LysC', 'LysN', 'GluC', 'ArgC'
                                                   'AspN',]),
              default=DEFAULTS['enzyme'], show_default=True,
              help="Enzyme used for digestion. Ignored with autorun.")
@click.option('-i', '--interval', type=int, default=DEFAULTS['interval'],
              help='''(Autorun only) Interval in seconds to wait between automatic
              checks for new files to group. Default is 1 hour.''')
@click.option('--ion-score', type=float, default=DEFAULTS['ion_score'], show_default=True,
              help='Ion score cutoff for a psm.')
@click.option('--ion-score-bins', nargs=3, type=float, default=DEFAULTS['ion_score_bins'], show_default=True,
              help='The three (ascending) IonScore cutoffs used to place PSMs in quality bins.')
@click.option('-l', '--labeltype', type=click.Choice(['none', 'SILAC', 'iTRAQ', 'TMT']),
              default=DEFAULTS['labeltype'], show_default=True, help='Type of label for this experiment.')
@click.option('-m', '--max-files', type=int, default=DEFAULTS['max_files'],
              help='(Autorun only) Maximum number of experiments to queue for autorun')
@click.option('--modi', type=int, default=DEFAULTS['modi'], show_default=True,
              help='Maximum modifications to allow on one peptide')
@click.option('-n', '--name', type=str, default=getuser(), show_default=True,
              help='Name associated with the search')
@click.option('-ntr', '--no-taxa-redistrib', is_flag=True, show_default=True,
              help='Disable redistribution based on individual taxons')
@click.option('-o', '--outdir', type=click.Path(file_okay=False), default=None,
              help='Output directory for files.')
@click.option('-p', '--psms-file', type=click.Path(exists=True, dir_okay=False),
              multiple=True,
              help='Tab deliminated file of psms to be grouped')
@click.option('--pipeline', type=click.Choice(['PD', 'MQ']), show_default=True,
              help='Pipeline used for generating PSM file.')
@click.option('--idg', type=int, default=DEFAULTS['idg'], show_default=True,
              help='PSM IDG cutoff value.')
@click.option('--pep', type=float, default=DEFAULTS['pep'], show_default=True,
              help='Posterior error probability cutoff')
@click.option('--qvalue', type=float, default=DEFAULTS['qvalue'], show_default=True,
              help='Cutoff q-value for a given psm.')
@click.option('--quant_source', type=click.Choice(['AUC', 'Intensity']), default=DEFAULTS['quant_source'],
              show_default=True, help='Cutoff q-value for a given psm.')
@click.option('-r', '--rawfiledir', type=click.Path(file_okay=False), default=DEFAULTS['rawfiledir'],
              show_default=True,
              help='''Directory to look for corresponding rawfiles for extra summary information.
              Note the raw files are not required for successful analysis.''')
@click.option('-s', '--configfile', type=click.Path(exists=True, dir_okay=False),
              help=''''Points to a specific configfile for use in the analysis.
              Note will automatically look for a `pygrouper_config.ini` in present directory if not specified''')
@click.option('-t', '--taxonid', type=int,
              help='Taxon ID associated with the database file')
@click.option('--zmin', type=int, default=DEFAULTS['zmin'], show_default=True,
              help='Minimum charge')
@click.option('--zmax', type=int, default=DEFAULTS['zmax'], show_default=True,
              help='Maximum charge')
def run(autorun, contaminants, database, enzyme, interval, ion_score, ion_score_bins, labeltype, max_files,
        modi, name, no_taxa_redistrib, outdir, psms_file, pipeline, idg, pep, qvalue, quant_source,
        rawfiledir, configfile, taxonid, zmin, zmax):
    """Run PyGrouper"""
    if not all([database, psms_file]) and not autorun:
        click.echo('No database or psms file entered, showing help and exiting...')
        click.echo(click.get_current_context().get_help())
        sys.exit(0)

    # sys.exit(0)
    config = parse_configfile(configfile)  # will parse if config file is specified or pygrouper_config.ini exists in PD
    if config is None:
        config = Config(name)
    INPUT_DIR = config.inputdir
    OUTPUT_DIR = outdir or config.outputdir or '.'
    RAWFILE_DIR = rawfiledir or config.rawfiledir
    LABELS = config.labels
    refseqs = config.refseqs
    filtervalues = config.filtervalues
    if filtervalues:
        click.echo('Using predefined config file')
    column_aliases = config.column_aliases
    gid_ignore_file = contaminants or config.contaminants
    if autorun:
        auto_grouper.interval_check(interval, INPUT_DIR, OUTPUT_DIR,
                                    max_files, rawfilepath=RAWFILE_DIR,
                                    refs=refseqs,
                                    column_aliases=column_aliases,
                                    gid_ignore_file=gid_ignore_file,
                                    labels=LABELS,)
    else:
        usrdatas = list()
        if not taxonid:
            taxonid = click.prompt('Enter taxon id', default=9606, type=int,)
        for ix, psmfile in enumerate(psms_file):
            usrdata = UserData(datafile=psmfile, indir=INPUT_DIR, outdir=OUTPUT_DIR, rawfiledir=RAWFILE_DIR,
                               no_taxa_redistrib=no_taxa_redistrib, labeltype=labeltype, addedby=name,
                               searchdb=database)
            try:
                rec, run, search = find_rec_run_search(psmfile)
                usrdata.recno, usrdata.runno, usrdata.searchno = int(rec), int(run), int(search)
            except AttributeError:  # regex search failed, just use a default
                usrdata.recno = ix+1  # default recno starts at 1
            usrdata.taxonid = taxonid
            refseqs[taxonid] = database
            INPUT_DIR, usrfile = os.path.split(Path(psmfile).resolve().__str__())
            usrdata.indir, usrdata.datafile = INPUT_DIR, usrfile
            usrdata.outdir = OUTPUT_DIR or Path(OUTPUT_DIR).resolve().__str__()
            # later on expected that datafile is separated from path
            usrdata.quant_source = quant_source
            usrdata.pipeline = pipeline
            if filtervalues: # if defined earlier from passed config file
                usrdata.filtervalues = filtervalues
                params = click.get_current_context().params
                for param in ('ion_score', 'ion_score_bins', 'qvalue', 'pep', 'idg', 'zmin', 'zmax', 'modi'):  # need to check for explictly passed options
                    if params[param] != DEFAULTS[param]:  # Explicitly over-write
                        usrdata.filtervalues[param] = params[param]
            else:
                usrdata.filtervalues['ion_score']      = ion_score
                usrdata.filtervalues['ion_score_bins'] = ion_score_bins
                usrdata.filtervalues['qvalue']         = qvalue
                usrdata.filtervalues['pep']            = pep
                usrdata.filtervalues['idg']            = idg
                usrdata.filtervalues['zmin']           = zmin
                usrdata.filtervalues['zmax']           = zmax
                usrdata.filtervalues['modi']           = modi



            usrdatas.append(usrdata)
        pygrouper.main(usrdatas=usrdatas,
                       inputdir=INPUT_DIR, outputdir=OUTPUT_DIR,
                       refs=refseqs, column_aliases=column_aliases,
                       gid_ignore_file=contaminants, labels=LABELS)

def find_rec_run_search(target):
    "Try to get record, run, and search numbers with regex of a target string with pattern \d+_\d+_\d+"
    _, target = os.path.split(target)  # ensure just searching on filename
    rec_run_search = re.compile(r'^\d+_\d+_\d+_')
    match = rec_run_search.search(target).group()
    recno = re.search(r'^\d+', match).group()
    recno_pat = re.compile('(?<={}_)\d+'.format(recno))
    runno = re.search(recno_pat, match).group()
    runno_pat = re.compile('(?<={}_{}_)\d+'.format(recno, runno))
    searchno = re.search(runno_pat, match).group()
    return recno, runno, searchno

if __name__ == '__main__':
    print('hi')
    # from click.testing import CliRunner
    # test()
