import os
import sys
import re
import unittest
import copy
import traceback
try:
    from unittest import mock
except ImportError:
    import mock
from io import StringIO
from click.testing import CliRunner
import numpy as np
import pandas as pd

from gpgrouper import gpgrouper
# from pygrouper import auto_grouper

from gpgrouper import cli
# from .. import cli
from gpgrouper.containers import UserData
# from ..containers import UserData
from RefProtDB.utils import _fasta_dict_from_file, fasta_dict_from_file


stout = StringIO()  # capture all of the click.echos here

# INPUT_DIR, PSMS_FILE = os.path.split('./testdata/two_uniques/test_input_all_02.tab')
BASEDIR, _ = os.path.split(os.path.abspath(__file__))

INPUT_DIR = os.path.join(BASEDIR, 'testdata/two_uniques')
PSMS_FILE = os.path.join(BASEDIR, 'testdata/two_uniques/test_input_all_02.tab')
OUTPUT_DIR = RAWFILE_DIR = INPUT_DIR
# REFSEQ_FILE = os.path.join(BASEDIR, 'testdata/two_uniques/refseq_02.tab')
REFSEQ_FILE = os.path.join(BASEDIR, 'testdata/two_uniques/refseq_02.fa')
CONFIG_FILE = os.path.join(BASEDIR, '../pygrouper_config.ini')

# REQUIRED_HEADERS = pygrouper.REQUIRED_HEADERS


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
        for dir_ in ('.', './testdata/', './testdata/two_uniques/'):
            for f in os.listdir(dir_):
                if f.startswith('10101_1_1') or f.startswith('1_1_1'):
                    os.remove(os.path.join(dir_, f))



    @mock.patch('gpgrouper.pygrouper.main')
    @mock.patch('sys.stdout', stdout)
    def test_call_run(self, main):
        self.longMessage = True
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run'])
        self.assertEqual(response.exit_code, 1, msg=response.output)

    # @mock.patch('sys.stdout', devnull)
    @mock.patch('gpgrouper.pygrouper.main')
    def test_call_run(self, main):
        self.longMessage = True
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
                                           '--psms-file', PSMS_FILE])
        self.assertEqual(response.exit_code, 0, msg=response.output)

    @mock.patch('gpgrouper.gpgrouper.main')
    def test_simple_input(self, main):
        """Test calling of `gpgrouper.main` from `gpgrouper run` with valid database, psms file, and taxonid"""
        self.longMessage = True
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
                                           '--psms-file', PSMS_FILE,
                                           '--outdir', INPUT_DIR,
                                           '--configfile', CONFIG_FILE,
                                           '--taxonid', 9606])
        self.assertTrue(main.called, msg=response.output)

    # @mock.patch('pygrouper.pygrouper.main')
    # def test_two_inputs(self, main):
    #     """Test if we can run with two inputs files"""
    #     runner = CliRunner()
    #     response = runner.invoke(cli.cli, ['run', '--database', REFSEQ_FILE,
    #                                        '--psms-file', PSMS_FILE,
    #                                        '--psms-file', PSMS_FILE,
    #                                        '--configfile', CONFIG_FILE,
    #                                        '--taxonid', 9606,
    #     ])
    #     self.assertEqual(len(main.call_args[-1]['usrdatas']), 2,
    #                      msg=response.exception)

    def test_load_fasta(self):
        with open(REFSEQ_FILE, 'r') as f:
            total_records = f.read().count('>')
        out = gpgrouper.load_fasta(REFSEQ_FILE)
        self.assertEqual(total_records, len(out))

    def test_invalid_fasta(self):
        orig = gpgrouper.fasta_dict_from_file
        gpgrouper.__dict__['fasta_dict_from_file'] = _fasta_dict_from_file
        fasta = StringIO()
        fasta.write(u'>gi|1234\nXXXXXX')
        fasta.seek(0)
        try:
            with self.assertRaises(ValueError):
                gpgrouper.load_fasta(fasta)
        finally:
            gpgrouper.__dict__['fasta_dict_from_file'] = orig



def get_sample_data():
        runner = CliRunner()
        with mock.patch('gpgrouper.gpgrouper.main') as main:
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
        for dir_ in ('.', './testdata/', './testdata/two_uniques/'):
            for f in os.listdir(dir_):
                if f.startswith('10101_1_1') or f.startswith('1_1_1'):
                    os.remove(os.path.join(dir_, f))


    def test_called_match(self):
        "Test if match function gets called (stage_match is working properly)"
        sample = get_sample_data()
        with mock.patch('gpgrouper.gpgrouper._match') as _match:
            gpgrouper.match(sample['usrdatas'], sample['refs'])
            self.assertTrue(_match.called)

    def test_called_matched_multiple(self):
        """Test if match function gets called twice with multiple taxonids
        (stage_match is working properly)"""
        sample = get_sample_data()
        usrdata1 = copy.copy(sample['usrdatas'][-1])
        usrdata1.taxonid = 10090
        sample['refs'][10090] = sample['refs'][9606]
        sample['usrdatas'].append(usrdata1)
        with mock.patch('gpgrouper.gpgrouper._match') as _match:
            gpgrouper.match(sample['usrdatas'], sample['refs'])
            self.assertEqual(_match.call_count, 2)

    def test_matcher(self):
        """Test functionality of matcher,
        does not validate any data, just if it dones't fail """
        sample = get_sample_data()
        usrdatas = sample['usrdatas']
        gpgrouper.set_up(usrdatas, sample['column_aliases'])
        _ = gpgrouper.match(usrdatas, sample['refs'])
        usrdata = usrdatas[-1]
        outcols = ['GeneIDs_All', 'GeneIDCount_All', 'TaxonIDs_All', 'TaxonIDCount_All',
                   'ProteinGIs_All', 'ProteinGICount_All']
        for outcol in outcols:
            self.assertIn(outcol, usrdata.df.columns,
                          msg='Matcher not returning correct columns')

_logfile = re.compile('gpGrouper_v.*.log')

class TestMin(unittest.TestCase):
    """Test with the minimum required data"""

    _dirname = os.path.dirname(os.path.realpath(__file__))

    FASTA = os.path.join(_dirname, './testdata/minfasta.fa')
    PSMS  = os.path.join(_dirname, './testdata/minfile.tab')

    stdout = sys.stdout
    stderr = sys.stderr

    def setUp(self):
        sys.stdout = StringIO()
        sys.stderr = StringIO()
        # fasta_h = '\t'.join(['TaxonID', 'HomologeneID', 'GeneID', 'ProteinGI', 'FASTA'])
        # fasta_values = '\t'.join(['9606', '', '1234', '12345678', 'AAAAAAA'])
        # fasta_file = '\n'.join([fasta_h, fasta_values])

        fasta_h = '>geneid|{gid}|ref|{ref}|taxon|{taxon}|gi|{gi}|homologene|{hid}| {desc}'
        fasta_values = ''

        r = dict(gid='1234', ref='NP_XX', taxon='9606', gi='12345678', hid='', desc=' \n')
        fasta_values += fasta_h.format(**r)
        fasta_values += 'AAAAAAA\n'
        fasta_file = fasta_values

        headers = '\t'.join(['Sequence', 'Modifications', 'PrecursorArea',
                             'Charge', 'IonScore', 'q_value', 'PEP', 'SpectrumFile',
                             'RTmin', 'DeltaMassPPM'])
        values = '\t'.join(['AAAAAAA', ' ', '12345', '2', '50', '0.00', '0.00', 'file1.raw', '10', '1'])
        small_file = '\n'.join([headers, values])

        with open(self.FASTA, 'w') as fasta, open(self.PSMS, 'w') as file_:
            fasta.write(fasta_file)
            file_.write(small_file)

    def tearDown(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        os.remove(self.FASTA)
        os.remove(self.PSMS)
        for dir_ in ('.', './testdata/', './testdata/two_uniques/'):
            for f in os.listdir(dir_):
                if f.startswith('10101_1_1') or f.startswith('1_1_1'):
                    os.remove(os.path.join(dir_, f))

    # def test_invalid_file(self):
    #     runner = CliRunner()
    #     df = (pd.read_table(self.PSMS)
    #           .drop(['PrecursorArea'], axis=1))

    #     with self.assertRaises(ValueError):
    #         pygrouper.check_required_headers(df)

    def test_minimum_cols(self):
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', self.FASTA,
                                           '--psms-file', self.PSMS,
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

    FASTA = os.path.join(_dirname, './testdata/tmt_fa.tab')
    PSMS  = os.path.join(_dirname, './testdata/10101_1_1_tmt_fa.tab')

    stdout = sys.stdout
    stderr = sys.stderr

    set2s = [2345, 1234, 678, 567]
    ratios = {'9606' : 4/5.0,
              '10090' : 1/5.0}

    def setUp(self):
        sys.stdout = StringIO()
        sys.stderr = StringIO()

        # fasta_h = '\t'.join(['TaxonID', 'HomologeneID', 'GeneID', 'ProteinGI', 'FASTA'])
        # fasta_values = '\t'.join(['9606', '', '1234', '12345678', 'AAAAAAAKBBBBBBBKCCCCCCCK\n'])
        # fasta_values += '\t'.join(['9606', '', '2345', '12345678', 'AAAAAAAKBBBBBBBKCCCCCCCK\n'])
        # fasta_values += '\t'.join(['9606', '', '111', '11111', 'CBACBACBAK\n'])
        # fasta_values += '\t'.join(['9606', '', '666', '46123', 'BBBBBBBKCCCCCCCK\n'])
        # fasta_values += '\t'.join(['10090', '', '567', '865431', 'AAAAAAAKBBBBBBBKCCCCCCCK\n'])
        # fasta_values += '\t'.join(['10090', '', '678', '865431', 'AAAAAAAKBBBBBBBKCCCCCCCK\n'])
        # fasta_values += '\t'.join(['10090', '', '999', '99999', 'ABCABCABCK\n'])
        # fasta_file = '\n'.join([fasta_h, fasta_values])

        fasta_h = '>geneid|{gid}|ref|{ref}|taxon|{taxon}|gi|{gi}|homologene|{hid}| {desc}'
        fasta_values = ''

        r = dict(gid='1234', ref='NP_XX', taxon='9606', gi='12345678', hid='1', desc=' \n')
        fasta_values += fasta_h.format(**r)
        fasta_values += 'AAAAAAAKBBBBBBBKCCCCCCCK\n'

        r = dict(gid='2345', ref='NP_XX', taxon='9606', gi='12345678', hid='1', desc=' \n')
        fasta_values += fasta_h.format(**r)
        fasta_values += 'AAAAAAAKBBBBBBBKCCCCCCCK\n'

        r = dict(gid='111', ref='NP_XX', taxon='9606', gi='11111', hid='1', desc=' \n')
        fasta_values += fasta_h.format(**r)
        fasta_values += 'CBACBACBAK\n'

        r = dict(gid='666', ref='NP_XX', taxon='9606', gi='46123', hid='1', desc=' \n')
        fasta_values += fasta_h.format(**r)
        fasta_values += 'BBBBBBBKCCCCCCCK\n'

        r = dict(gid='567', ref='NP_XX', taxon='10090', gi='865431', hid='1', desc=' \n')
        fasta_values += fasta_h.format(**r)
        fasta_values += 'AAAAAAAKBBBBBBBKCCCCCCCK\n'

        r = dict(gid='678', ref='NP_XX', taxon='10090', gi='865431', hid='1', desc=' \n')
        fasta_values += fasta_h.format(**r)
        fasta_values += 'AAAAAAAKBBBBBBBKCCCCCCCK\n'

        r = dict(gid='999', ref='NP_XX', taxon='10090', gi='99999', hid='1', desc=' \n')
        fasta_values += fasta_h.format(**r)
        fasta_values += 'ABCABCABCK\n'

        fasta_file = fasta_values


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
                                   'cBACBACBAK', 'aBCaBCABCK'],
                 Modifications  = ['N-Term(TMT6plex)'] * n,
                 PrecursorArea  = [np.multiply(x, np.power(10, 2))
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
        d['Modifications'][-1] += '; A4(Phospho)'

        df = pd.DataFrame(d)
        df.to_csv(self.PSMS, sep='\t', index=False)

        with open(self.FASTA, 'w') as fasta:
            fasta.write(fasta_file)

    def tearDown(self):
        # return
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        os.remove(self.FASTA)
        os.remove(self.PSMS)
        for dir_ in ('.', './testdata/', './testdata/two_uniques/'):
            for f in os.listdir(dir_):
                if f.startswith('10101_1_1') or f.startswith('1_1_1'):
                    os.remove(os.path.join(dir_, f))


        for f in os.listdir('.') + os.listdir('./testdata/') + os.listdir('./testdata/two_uniques'):
            if f.startswith('10101_1_1') or f.startswith('1_1_1'):
                try:
                    os.remove(f)
                except (FileNotFoundError):
                    pass

    def test_phospho(self):
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', self.FASTA,
                                           '--psms-file', self.PSMS,
                                           '--taxonid', 9606,
                                           '--outdir', './testdata',
                                           '--configfile', CONFIG_FILE,
                                           '--labeltype', 'none',
                                           '--phospho'
                                           ])

        # print full traceback if there is a failure
        self.assertEqual(0, response.exit_code,
                         msg='\n{}\n{!r}'.format(''.join(traceback.format_tb(response.exc_info[-1])),
                                                                             response.exc_info[1])
                         )

        df = pd.read_table('./testdata/10101_1_1_none_e2g.tab').fillna(0)
        self.assertGreater( df.query('GeneID==999').AreaSum_dstrAdj.values[0], 0)

        for gid in df.GeneID.unique():  # rest should be zero
            if gid == 999:
                continue
            self.assertEqual(df.loc[df.GeneID==gid, 'AreaSum_dstrAdj'].values[0],
                               0)

            self.assertGreater(df.loc[df.GeneID==gid, 'PeptideCount'].values[0],
                               0)


    def test_labelfree(self):
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', self.FASTA,
                                           '--psms-file', self.PSMS,
                                           '--taxonid', 9606,
                                           '--outdir', './testdata',
                                           '--configfile', CONFIG_FILE,
                                           '--labeltype', 'none',
                                           ])
        logf = open('./testdata/10101_1_1_none.log').read()
        # print full traceback if there is a failure
        self.assertEqual(0, response.exit_code,
                         msg='\n{}\n{!r}'.format(''.join(traceback.format_tb(response.exc_info[-1])),
                                                                             response.exc_info[1])
                         )
        dstrAdj = 'AreaSum_dstrAdj'
        maxarea = 'AreaSum_max'
        df = pd.read_table('./testdata/10101_1_1_none_e2g.tab')
        self.assertEqual(True, all(df[ df.GeneID.isin(self.set2s)]['IDSet'] == 2))
        self.assertEqual(True, all(df.query('IDSet==3')['AreaSum_dstrAdj']==0))
        for tid in 9606, 10090:
            q = 'IDSet==2 & TaxonID == {}'.format(tid)
            subdf = df.query(q)
            tot = len(subdf)

            # all(subdf[dstrAdj]==(subdf[maxarea]*self.ratios[str(tid)])/tot),
            msg = '\n'+str(subdf[['GeneID', 'TaxonID', 'IDSet', dstrAdj, maxarea]]) +'\n' +\
                  str(subdf[maxarea]*self.ratios[str(tid)]/tot) +\
                  '\n'+ str(self.ratios[str(tid)]) +\
                  '\n'+ logf
                  # '\n' + str(df) +\
            # self.assertEqual(True,
            #                  all(subdf[dstrAdj]==(subdf[maxarea]*self.ratios[str(tid)])/tot),
            #                  msg=msg
            #                  # msg=self.ratios[str(tid)]
            # )
            np.testing.assert_almost_equal(
                subdf[dstrAdj],
                (subdf[maxarea]*self.ratios[str(tid)]/tot).values,
                decimal=7,
                err_msg=msg,
                             # msg='\n'+str(subdf[[dstrAdj, maxarea]]) +'\n'+ str(self.ratios[str(tid)])
                             # msg=subdf
                             # msg=self.ratios[str(tid)]
            )




# @mock.patch('sys.stdout', devnull)
class TestFull(unittest.TestCase):
    """Test a full runthrough to make sure we can go from start to finish"""

    stdout = sys.stdout
    stderr = sys.stderr

    FASTA = REFSEQ_FILE
    PSMS = PSMS_FILE

    def setUp(self):
        sys.stdout = StringIO()
        sys.stderr = StringIO()
        # sys.stdout = self.stdout
        # sys.stderr = self.stderr


    def tearDown(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        for dir_ in ('.', './testdata/', './testdata/two_uniques/'):
            for f in os.listdir(dir_):
                if f.startswith('10101_1_1') or f.startswith('1_1_1'):
                    os.remove(os.path.join(dir_, f))

    def test_runthrough(self):
        # self.longMessage = True
        runner = CliRunner()
        response = runner.invoke(cli.cli, ['run', '--database', self.FASTA,
                                           '--psms-file', self.PSMS,
                                           '--taxonid', 9606,
                                           '--outdir', INPUT_DIR,
                                           '--configfile', CONFIG_FILE,
                                           ])
        # print full traceback if there is a failure
        self.assertEqual(0, response.exit_code,
                         msg='\n{}\n{!r}'.format(''.join(traceback.format_tb(response.exc_info[-1])),
                                                                             response.exc_info[1])
                         )
    def test_proper_columns_psms(self):

        data_cols = gpgrouper.DATA_COLS
        try:
            data_cols.remove('LastScan')
        except ValueError:
            pass

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
            self.assertTrue(col in output.columns, msg='{} not found in data file'.format(col))


    # @mock.patch('pygrouper.auto_grouper._update_database')
    # def test_update_db(self, func):
    #     """Test if update database is successfully called"""
    #     sample_data = get_sample_data()
    #     # print(sample_data)
    #     auto_grouper.run_and_update(**sample_data)
    #     self.assertEqual(func.call_count, 1)



if __name__ == '__main__':
    unittest.main()
    for f in os.listdir('.') + os.listdir('./testdata/') + os.listdir('./testdata/two_uniques/'):
        print(f)
        if f.startswith('10101_1_1') or f.startswith('1_1_1'):
            # try:
            os.remove(f)
            # except (PermissionError, FileNotFoundError):
            #     pass
