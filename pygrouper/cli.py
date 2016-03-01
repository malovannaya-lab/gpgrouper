import os
import click
from pygrouper.manual_tests import test as manual_test

__author__ = 'Alexander B. Saltzman'
__copyright__ = 'Copyright January 2016'
__credits__ = ['Alexander B. Saltzman', 'Anna Malovannaya']
__license__ = 'MIT'
__version__ = '0.1.015'
__maintainer__ = 'Alexander B. Saltzman'
__email__ = 'saltzman@bcm.edu'


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
class Config(object):

    def __init__(self):
        self.email = None
        self.outfile = '-'

pass_config = click.make_pass_decorator(Config, ensure=True)

@click.group()
def cli():
    pass

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-q', '--quick', is_flag=True, help='Run a single small file')
def test(quick, ):
    """Test pygrouper with some pre-existing data."""
    manual_test.runtest(quick)
