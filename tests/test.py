import os
import sys
import re
import unittest
import copy
import traceback
from unittest import mock
from io import StringIO
from click.testing import CliRunner
import numpy as np
import pandas as pd

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
# devnull = open(os.devnull, 'w')

class InputTest(unittest.TestCase):
    """Test ability to call `pygrouper run`"""

    stdout = sys.stdout
    stderr = sys.stderr

    def setUp(self):
        sys.stdout = StringIO()
        sys.stderr = StringIO()

    def tearDown(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr



    # @mock.patch('sys.stdout', devnull)
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
                                           '--configfile', CONFIG_FILE,
                                           '--taxonid', 9606])
        self.assertTrue(main.called, msg=response.output)

    @mock.patch('pygrouper.pygrouper.main')
    def test_two_inputs(self, main):
        """Test if we can run with two inputs files"""
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
                                           '--psms-file', PSMS_FILE,
                                           '--psms-file', PSMS_FILE,
                                           '--configfile', CONFIG_FILE,
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
                                               '--outdir', INPUT_DIR,
                                               '--configfile', CONFIG_FILE,
            ]
                                     )
            # response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
            #                                    '--psms-file', PSMS_FILE,
            #                                    '--taxonid', 10090]
            #                          )
        return main.call_args[-1]

# @mock.patch('sys.stdout', devnull)
class MatchTest(unittest.TestCase):


    stdout = sys.stdout
    stderr = sys.stderr

    def setUp(self):
        sys.stdout = StringIO()
        sys.stderr = StringIO()

    def tearDown(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr

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
        usrdatas = pygrouper.set_up(sample['usrdatas'], sample['column_aliases'])
        usrdatas, _ = pygrouper.match(usrdatas, sample['refs'])
        usrdata = usrdatas[-1]
        outcols = ['GeneList', 'GeneCount', 'TaxonIDList', 'TaxonCount',
                   'ProteinList', 'ProteinCount']
        for outcol in outcols:
            self.assertIn('psm_' + outcol, usrdata.df.columns,
                          msg='Matcher not returning correct columns')

_logfile = re.compile('Pygrouper_v.*.log')

class TestMin(unittest.TestCase):
    """Test with the minimum required data"""

    _dirname = os.path.dirname(os.path.realpath(__file__))

    MIN_FASTA = os.path.join(_dirname, './testdata/minfasta.tab')
    MIN_FILE  = os.path.join(_dirname, './testdata/minfile.tab')

    stdout = sys.stdout
    stderr = sys.stderr

    def setUp(self):
        sys.stdout = StringIO()
        sys.stderr = StringIO()
        fasta_h = '\t'.join(['TaxonID', 'HomologeneID', 'GeneID', 'ProteinGI', 'FASTA'])
        fasta_values = '\t'.join(['9606', '', '1234', '12345678', 'AAAAAAA'])
        fasta_file = '\n'.join([fasta_h, fasta_values])

        headers = '\t'.join(['Sequence', 'Modifications', 'PrecursorArea',
                             'Charge', 'IonScore', 'q_value', 'PEP', 'SpectrumFile',
                             'RTmin', 'DeltaMassPPM'])
        values = '\t'.join(['AAAAAAA', ' ', '12345', '2', '50', '0.00', '0.00', 'file1.raw', '10', '1'])
        small_file = '\n'.join([headers, values])

        with open(self.MIN_FASTA, 'w') as fasta, open(self.MIN_FILE, 'w') as file_:
            fasta.write(fasta_file)
            file_.write(small_file)

    def tearDown(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        os.remove(self.MIN_FASTA)
        os.remove(self.MIN_FILE)
        for f in os.listdir('.'):
            if f.startswith('1_1_1'):
                os.remove(f)

    def test_minimum_cols(self):
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', self.MIN_FASTA,
                                           '--psms-file', self.MIN_FILE,
                                           '--taxonid', 9606,
                                           '--outdir', '.',
                                           '--configfile', CONFIG_FILE,
                                           ])
        # print full traceback if there is a failure
        self.assertEqual(0, response.exit_code,
                         msg='\n{}\n{!r}'.format(''.join(traceback.format_tb(response.exc_info[-1])),
                                                                             response.exc_info[1])
                         )

class TestAreaTMT(unittest.TestCase):
    """ Testing for dealing with Areas for TMT and Labelfree data"""
    _dirname = os.path.dirname(os.path.realpath(__file__))

    TMT_FASTA = os.path.join(_dirname, './testdata/tmt_fa.tab')
    TMT_FILE  = os.path.join(_dirname, './testdata/10101_1_1_tmt_fa.tab')

    stdout = sys.stdout
    stderr = sys.stderr

    set2s = [2345, 1234, 678, 567]
    ratios = {'9606' : 4/5.0,
              '10090' : 1/5.0}

    def setUp(self):
        sys.stdout = StringIO()
        sys.stderr = StringIO()

        fasta_h = '\t'.join(['TaxonID', 'HomologeneID', 'GeneID', 'ProteinGI', 'FASTA'])
        fasta_values = '\t'.join(['9606', '', '1234', '12345678', 'AAAAAAAKBBBBBBBKCCCCCCCK\n'])
        fasta_values += '\t'.join(['9606', '', '2345', '12345678', 'AAAAAAAKBBBBBBBKCCCCCCCK\n'])
        fasta_values += '\t'.join(['9606', '', '111', '11111', 'CBACBACBAK\n'])
        fasta_values += '\t'.join(['9606', '', '666', '46123', 'BBBBBBBKCCCCCCCK\n'])

        fasta_values += '\t'.join(['10090', '', '567', '865431', 'AAAAAAAKBBBBBBBKCCCCCCCK\n'])
        fasta_values += '\t'.join(['10090', '', '678', '865431', 'AAAAAAAKBBBBBBBKCCCCCCCK\n'])
        fasta_values += '\t'.join(['10090', '', '999', '99999', 'ABCABCABCK\n'])
        fasta_file = '\n'.join([fasta_h, fasta_values])


        n = 3 # number of PSMs
        # d = dict(Sequence       = ['aAAAAAAK', 'bBBBBBBK', 'bBBBBBBK'],
        #          Modifications  = ['N-Term(TMT6plex)'] * n,
        #          PrecursorArea  = [100] * n,
        #          Charge         = [2] * n,
        #          IonScore       = [50, 50, 30],
        #          q_value        = [0] * n,
        #          PEP            = [0] * n,
        #          SpectrumFile   = ['File1.raw'] * n,
        #          RTmin          = [1, 2, 2.2],
        #          DeltaMassPPM   = [0] * n,
        #          TMT_126        = [2, 2, 4],
        #          TMT_131        = [4, 4, 2]
        # )
        n = 5

        d = dict(Sequence       = ['aAAAAAAK', 'bBBBBBBK', 'bBBBBBBK',
                                   'cBACBACBAK', 'aBCABCABCK'],
                 Modifications  = ['N-Term(TMT6plex)'] * n,
                 PrecursorArea  = [np.multiply(x, np.power(10, 8))
                                   for x in [100, 100, 100, 80, 20]],
                 Charge         = [2] * n,
                 IonScore       = [50]* n,
                 q_value        = [0] * n,
                 PEP            = [0] * n,
                 SpectrumFile   = ['File1.raw'] * n,
                 RTmin          = [1.1, 1, 1, 2.2, 2.3],
                 DeltaMassPPM   = [0] * n,
                 TMT_126        = [2, 2, 2, 4, 4],
                 TMT_131        = [4, 4, 4, 4, 4]
        )

        df = pd.DataFrame(d)
        df.to_csv(self.TMT_FILE, sep='\t', index=False)

        with open(self.TMT_FASTA, 'w') as fasta:
            fasta.write(fasta_file)

    def tearDown(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        os.remove(self.TMT_FASTA)
        os.remove(self.TMT_FILE)
        for f in os.listdir('.'):
            if f.startswith('10101_1_1'):
                os.remove(f)
        for f in os.listdir('./testdata'):
            if f.startswith('10101_1_1'):
                os.remove(os.path.join('./testdata', f))

    def test_labelfree(self):
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', self.TMT_FASTA,
                                           '--psms-file', self.TMT_FILE,
                                           '--taxonid', 9606,
                                           '--outdir', './testdata',
                                           '--configfile', CONFIG_FILE,
                                           '--labeltype', 'none',
                                           ])

        # print full traceback if there is a failure
        self.assertEqual(0, response.exit_code,
                         msg='\n{}\n{!r}'.format(''.join(traceback.format_tb(response.exc_info[-1])),
                                                                             response.exc_info[1])
                         )
        dstrAdj = 'e2g_nGPArea_Sum_dstrAdj'
        maxarea = 'e2g_nGPArea_Sum_max'
        df = pd.read_table('./testdata/10101_1_1_none_0_e2g.tab')
        print(df.columns)
        self.assertEqual(True, all(df[ df.e2g_GeneID.isin(self.set2s)]['e2g_IDSet'] == 2))
        self.assertEqual(True, all(df.query('e2g_IDSet==3')['e2g_nGPArea_Sum_dstrAdj']==0))
        for tid in 9606, 10090:
            q = 'e2g_IDSet==2 & e2g_TaxonID == {}'.format(tid)
            subdf = df.query(q)
            tot = len(subdf)
            self.assertEqual(True,
                             all(subdf[dstrAdj]==(subdf[maxarea]*self.ratios[str(tid)])/tot),
                             # msg=subdf[[dstrAdj, maxarea]]
                             msg=self.ratios[str(tid)]
            )




# @mock.patch('sys.stdout', devnull)
class TestFull(unittest.TestCase):
    """Test a full runthrough to make sure we can go from start to finish"""

    stdout = sys.stdout
    stderr = sys.stderr

    def setUp(self):
        sys.stdout = StringIO()
        sys.stderr = StringIO()


    def tearDown(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        for f in os.listdir('.'):
            if _logfile.match(f):
                os.remove(f)
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
                                           '--configfile', CONFIG_FILE,
                                           ])
        # print full traceback if there is a failure
        self.assertEqual(0, response.exit_code,
                         msg='\n{}\n{!r}'.format(''.join(traceback.format_tb(response.exc_info[-1])),
                                                                             response.exc_info[1])
                         )
    def test_proper_columns_e2g(self):
        data_cols = ['psm_EXPRecNo', 'psm_EXPRunNo', 'psm_EXPSearchNo',
                     'psm_EXPTechRepNo', 'Sequence',
                     'PSMAmbiguity', 'Modifications', 'ActivationType',
                     'DeltaScore', 'DeltaCn', 'Rank', 'SearchEngineRank',
                     'PrecursorArea', 'q_value', 'PEP',
                     'IonScore', 'MissedCleavages',
                     'IsolationInterference', 'IonInjectTime',
                     'Charge', 'mzDa', 'MHDa',
                     'DeltaMassDa', 'DeltaMassPPM', 'RTmin',
                     'FirstScan', 'MSOrder', 'MatchedIons',
                     'SpectrumFile', 'psm_AddedBy', 'psm_oriFLAG',
                     'psm_CreationTS', 'psm_ModificationTS', 'psm_GeneID',
                     'psm_GeneList', 'psm_GeneCount', 'psm_ProteinGI',
                     'psm_ProteinList', 'psm_ProteinCount',
                     'psm_HID', 'psm_HIDList', 'psm_HIDCount',
                     'psm_TaxonID', 'psm_TaxonIDList', 'psm_TaxonCount',
                     'psm_PSM_IDG', 'psm_SequenceModi',
                     'psm_SequenceModiCount', 'psm_LabelFLAG',
                     'psm_PeptRank', 'psm_AUC_UseFLAG', 'psm_PSM_UseFLAG',
                     'psm_Peak_UseFLAG', 'psm_SequenceArea', 'psm_PrecursorArea_dstrAdj']
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
                                           '--psms-file', PSMS_FILE,
                                           '--taxonid', 9606,
                                           '--outdir', INPUT_DIR,
                                           '--configfile', CONFIG_FILE,
                                           ])
        # output = pd.read_table(os.path.join(INPUT_DIR, '1_1_1_none_0_e2g.tab'))
        output = pd.read_table(os.path.join(INPUT_DIR, '1_1_1_none_psms.tab'))
        # print(output.columns)
        for col in data_cols:
            self.assertTrue(col in output.columns, msg='{} not found in e2g file'.format(col))


    @mock.patch('pygrouper.auto_grouper._update_database')
    def test_update_db(self, func):
        """Test if update database is successfully called"""
        sample_data = get_sample_data()
        # print(sample_data)
        auto_grouper.run_and_update(**sample_data)
        self.assertEqual(func.call_count, 1)



if __name__ == '__main__':
    unittest.main()
