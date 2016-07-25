import os
import sys
import unittest
import copy
import traceback
from unittest import mock
from io import StringIO
from click.testing import CliRunner

from pygrouper import pygrouper, auto_grouper
from pygrouper import cli
from pygrouper.containers import UserData


stout = StringIO()  # capture all of the click.echos here

# INPUT_DIR, PSMS_FILE = os.path.split('./testdata/two_uniques/test_input_all_02.tab')
BASEDIR, _ = os.path.split(os.path.abspath(__file__))

INPUT_DIR = os.path.join(BASEDIR, 'testdata/two_uniques')
PSMS_FILE = os.path.join(BASEDIR,'testdata/two_uniques/test_input_all_02.tab')
OUTPUT_DIR = RAWFILE_DIR = INPUT_DIR
REFSEQ_FILE = os.path.join(BASEDIR, 'testdata/two_uniques/refseq_02.tab')
CONFIG_FILE = os.path.join(BASEDIR, '../pygrouper_config.ini')

# usrdata = UserData(datafile=PSMS_FILE, indir=INPUT_DIR, outdir=OUTPUT_DIR, rawfiledir=RAWFILE_DIR,
                   # no_taxa_redistrib=0, labeltype='None', addedby='test', searchdb=REFSEQ_FILE)
devnull = open(os.devnull, 'w')

class InputTest(unittest.TestCase):
    """Test ability to call `pygrouper run`"""

    @mock.patch('sys.stdout', devnull)
    @mock.patch('pygrouper.pygrouper.main')
    def test_call_run(self, main):
        self.longMessage = True
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE])
        self.assertEqual(response.exit_code, 0, msg=response.output)

    @mock.patch('pygrouper.pygrouper.main')
    def test_simple_input(self, main):
        """Test calling of `pygrouper.main` from `pygrouper run` with valid database, psms file, and taxonid"""
        self.longMessage = True
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
                                           '--psms-file', PSMS_FILE,
                                           '--outdir', INPUT_DIR,
                                           '--taxonid', 9606])
        self.assertTrue(main.called, msg=response.output)

    @mock.patch('pygrouper.pygrouper.main')
    def test_two_inputs(self, main):
        """Test if we can run with two inputs files"""
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
                                           '--psms-file', PSMS_FILE,
                                           '--psms-file', PSMS_FILE,
                                           '--taxonid', 9606,
        ])
        self.assertEqual(len(main.call_args[-1]['usrdatas']), 2,
                         msg=response.exception)

def get_sample_data():
        runner = CliRunner()
        with mock.patch('pygrouper.pygrouper.main') as main:
            response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
                                               '--psms-file', PSMS_FILE,
                                               '--taxonid', 9606,
                                               '--configfile', CONFIG_FILE
            ]
                                     )
            # response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
            #                                    '--psms-file', PSMS_FILE,
            #                                    '--taxonid', 10090]
            #                          )
        return main.call_args[-1]

@mock.patch('sys.stdout', devnull)
class MatchTest(unittest.TestCase):

    def test_called_match(self):
        "Test if match function gets called (stage_match is working properly)"
        sample = get_sample_data()
        with mock.patch('pygrouper.pygrouper._match') as _match:
            pygrouper.match(sample['usrdatas'], sample['refs'])
            self.assertTrue(_match.called)

    def test_called_matched_multiple(self):
        """Test if match function gets called twice with multiple taxonids
        (stage_match is working properly)"""
        sample = get_sample_data()
        usrdata1 = copy.copy(sample['usrdatas'][-1])
        usrdata1.taxonid = 10090
        sample['refs'][10090] = sample['refs'][9606]
        sample['usrdatas'].append(usrdata1)
        with mock.patch('pygrouper.pygrouper._match') as _match:
            pygrouper.match(sample['usrdatas'], sample['refs'])
            self.assertEqual(_match.call_count, 2)

    def test_matcher(self):
        """Test functionality of matcher,
        does not validate any data, just if it dones't fail """
        sample = get_sample_data()
        usrdatas, _ = pygrouper.match(sample['usrdatas'], sample['refs'])
        usrdata = usrdatas[-1]
        outcols = ['GeneList', 'GeneCount', 'TaxonIDList', 'TaxonCount',
                   'ProteinList', 'ProteinCount']
        for outcol in outcols:
            self.assertIn('psm_' + outcol, usrdata.df.columns,
                          msg='Matcher not returning correct columns')


@mock.patch('sys.stdout', devnull)
class TestFull(unittest.TestCase):
    """Test a full runthrough to make sure we can go from start to finish"""

    def tearDown(self):
        for f in os.listdir(INPUT_DIR):
            if f.startswith('1_1_1'):
                os.remove(os.path.join(INPUT_DIR, f))

    def test_runthrough(self):
        # self.longMessage = True
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
                                           '--psms-file', PSMS_FILE,
                                           '--taxonid', 9606,
                                           '--outdir', INPUT_DIR,
                                           # '--configfile', CONFIG_FILE
                                           ])
        # print full traceback if there is a failure
        self.assertEqual(0, response.exit_code,
                         msg='\n{}\n{!r}'.format(''.join(traceback.format_tb(response.exc_info[-1])),
                                                                             response.exc_info[1])
                         )

    @mock.patch('pygrouper.auto_grouper._update_database')
    def test_update_db(self, func):
        """Test if update database is successfully called"""
        sample_data = get_sample_data()
        # print(sample_data)
        auto_grouper.run_and_update(**sample_data)
        self.assertEqual(func.call_count, 1)



if __name__ == '__main__':
    unittest.main()
