import sys
import os
import subprocess
from click.testing import CliRunner
from pygrouper import cli

from test import TestAreaTMT

BASEDIR, _ = os.path.split(os.path.abspath(__file__))
CONFIG_FILE = os.path.join(BASEDIR, '../pygrouper_config.ini')

class Interactive(TestAreaTMT):

    def setUp(self):
        super().setUp()
        sys.stdout = self.stdout
        sys.stderr = self.stderr


    def interact(self):
        self.setUp()
        call = ['pygrouper', 'run', '--database', self.TMT_FASTA,
                '--psms-file', self.TMT_FILE,
                '--taxonid', '9606',
                '--outdir', './testdata',
                '--configfile', CONFIG_FILE,
                '--labeltype', 'none',
        ]
        subprocess.call(call)


        runner = CliRunner()
        self.tearDown()

def main():
    Interactive().interact()


if __name__ == '__main__':
    main()
