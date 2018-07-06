import sys
import os
import subprocess
from click.testing import CliRunner
from pygrouper import cli

from test import TestAreaTMT, TestFull, TestMin

BASEDIR, _ = os.path.split(os.path.abspath(__file__))
CONFIG_FILE = os.path.join(BASEDIR, '../pygrouper_config.ini')

# class Interactive(TestAreaTMT):
class Interactive(TestFull):
# class Interactive(TestMin):

    def setUp(self):
        super().setUp()
        sys.stdout = self.stdout
        sys.stderr = self.stderr

    def tearDown(self):
        pass

    def interact(self):
        self.setUp()
        call = ['pygrouper', 'run', '--database', self.FASTA,
                '--psms-file', self.PSMS,
                '--taxonid', '9606',
                '--outdir', './testdata',
                '--configfile', CONFIG_FILE,
                # '--labeltype', 'TMT',
                '--labeltype', 'none',
        ]
        subprocess.call(call)


        runner = CliRunner()
        self.tearDown()

def main():
    Interactive().interact()


if __name__ == '__main__':
    main()
