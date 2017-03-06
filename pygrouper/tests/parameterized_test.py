"""Testing pygrouper with nose_parameterized
"""
import sys
import os
import re
from collections import namedtuple, defaultdict
from nose_parameterized import parameterized
import unittest
import pandas as pd
sys.path.append('..')
import pygrouper
pd.set_option('max_columns',99)


def get_sample_data():
    """Get all test data in all folters within './testdata' as a
    dictionary of directory paths -> list of files within directory
    """
    TestData = namedtuple('TestData', 'directory, inputdata, refseq, presplit, postsplit, genedata')
    sample_data = dict()
    for (dirpath, dirnames, files) in os.walk('./testdata'):
        if files: #not an empty list
            dirsplit = '\\' if '\\' in dirpath else '/'
            directory = dirpath.split(dirsplit)[1]
            inputdata = [x for x in files if 'input' in x]
            refseq    = [x for x in files if 'refseq' in x]
            presplit  = [x for x in files if 'presplit' in x]
            postsplit = [x for x in files if 'postsplit' in x]
            genedata  = [x for x in files if 'e2g' in x]
            if all(len(x) == 1 for x in [inputdata, refseq, presplit, postsplit, genedata]):
                sample_data[dirpath] = TestData._make([directory, inputdata[0], refseq[0],
                                                       presplit[0], postsplit[0], genedata[0]])

    return sample_data

def sortdf(df, sortby=None, **kwds):
    """Sort two dataframes equally before comparing.
    Note that indices may (likely) will not match.
    kwds are passed to pandas sort_values()"""
    if sortby is None:
        sortby = ['Sequence', 'FirstScan', 'IonScore', 'PEP']
    return df.sort_values(by=sortby, ascending=[True]*len(sortby), **kwds)

def update_df(df1, df2, reassignments=None, sortby=None):
    """Update input dataframe with previous output dataframe values to bring
    previous results along."""
    df1 = sortdf(df1, sortby=sortby)
    df2 = sortdf(df2, sortby=sortby)
    for r in reassignments:
        df1[r] = df2[r]
    df1.sort_index(inplace=True)
    df2.sort_index(inplace=True)
    print('updating df')
    return (df1, df2)

class GrouperTestCase(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        self.path        = None
        self.inputf      = None
        self.refseqf     = None
        self.presplitf   = None
        self.postsplitf  = None
        self.e2gfile     = None
        self.currentTest = None
        self.area_col    = None
        unittest.TestCase.__init__(self, methodName)

    def __str__(self):
        ''' Override this so that we know which instance it is '''
        return "%s(%s) (%s)" % (self._testMethodName, self.currentTest,
                                unittest.util.strclass(self.__class__))

    def seq_modiAssert(self):
        """Test sequence modi function"""
        usrdata = pd.DataFrame(self.presplit, copy=True)
        usrdata['Sequence'], usrdata['psm_SequenceModi'],\
            usrdata['psm_SequenceModiCount'],\
            usrdata['psm_LabelFLAG'] = \
                        list(zip(*usrdata.apply(lambda x :
                                pygrouper.seq_modi(x['Sequence'],
                                x['Modifications']),
                                axis=1)))
        sortby = ['FirstScan', 'IonScore', 'PEP', 'Sequence']
        # test Sequence
        for (ix1,calc), (ix2,result) in zip(sortdf(usrdata, sortby).Sequence.items(),
                                            sortdf(self.presplit, sortby).Sequence.items()):
            with self.subTest(calc=calc):
                self.assertEqual(calc, result,
                                 msg=('Error in sequence at input index {}, '
                                      'result index {}').format(ix1, ix2))
        # test SequenceModi
        for (ix1,calc), (ix2,result) in zip(sortdf(usrdata, sortby).psm_SequenceModi.items(),
                                            sortdf(self.presplit, sortby).psm_SequenceModi.items()):
            with self.subTest(calc=calc):
                self.assertEqual(calc, result,
                                 msg=('Error in sequence modi at input index {}, '
                                      'result index {}').format(ix1, ix2))
        # test SequenceModiCount
        for (ix1,calc), (ix2,result) in zip(sortdf(usrdata, sortby).psm_SequenceModiCount.items(),
                                            sortdf(self.presplit, sortby).psm_SequenceModiCount.items()):
            with self.subTest(calc=calc):
                self.assertEqual(calc, result,
                                 msg=('Error in modification count at input index {}, '
                                      'result index {}').format(ix1, ix2))
        # test LabelFLAG
        for (ix1,calc), (ix2,result) in zip(sortdf(usrdata, sortby).psm_LabelFLAG.items(),
                                            sortdf(self.presplit, sortby).psm_LabelFLAG.items()):
            with self.subTest(calc=calc):
                self.assertEqual(calc, result,
                                 msg=('Error in label flag at input index {}, '
                                      'result index {}').format(ix1, ix2))


    def gene_mapperAssert(self):
        """Test peptide to gene mapping.
        Lots of copy-pasted code as not currently defined in a stand-alone function.
        TODO: refactor this code in original pygrouper.py
        """

        # first re-assign previously calculated columns in test data to result data
        #reassign = ['Sequence', 'psm_SequenceModi', 'psm_SequenceModiCount', 'psm_LabelFLAG']
        #self.usrdata, self.presplit = update_df(self.usrdata, self.presplit, reassignments=reassign)
        usrdata = pd.DataFrame(self.presplit, copy=True)
        ref_reader = pygrouper.csv_reader(os.path.join(self.path, self.refseqf))
        prot = defaultdict(list)
        while True:
            try:
                row = next(ref_reader)

                fragments, fraglen = pygrouper.protease(row.fasta, minlen=7,
                                                        cutsites=['K', 'R'],
                                                        exceptions=['P'])
                for fragment in fragments:
                    prot[fragment].append(
                        pygrouper.RefseqInfo._make([row.taxonid, row.geneid,
                                                    row.homologeneid,row.proteingi,
                                                    fraglen]))
            except StopIteration:
                break
        usrdata['metadatainfo'] = ''
        usrdata['metadatainfo'] = usrdata.apply(lambda x:
                                                pygrouper.genematcher(x['Sequence'],
                                                            x['metadatainfo'],
                                                            prot), axis=1)
        # ==================== Populate gene info ================================ #
        gene_metadata = defaultdict(list)
        gene_taxon_dict = dict()
        for metadata in usrdata.metadatainfo:
            for data in metadata:
                gene_metadata[data.geneid].append((data.taxonid, data.homologeneid,
                                                   data.proteingi, data.genefraglen))
        for gene in gene_metadata:
            gene_taxon_dict[gene] = gene_metadata[gene][0][0]
        usrdata['psm_GeneList'], usrdata['psm_GeneCount'], \
            usrdata['psm_TaxonIDList'],\
            usrdata['psm_TaxonCount'], self.usrdata['psm_ProteinList'], \
            usrdata['psm_ProteinCount'] = list(zip(
                *usrdata.apply(lambda x : pygrouper.genelist_extractor(x['metadatainfo'],
                ),
                               axis=1)))
        # ==================== Populate gene info ================================ #

        def get_sequences(ix1, df1, ix2, df2):
            return (df1.loc[ix1]['Sequence'], df2.loc[ix2]['Sequence'])

        sortby = ['Sequence', 'FirstScan', 'IonScore', 'PEP']
        # test GeneList
        for (ix1,calc), (ix2,result) in zip(sortdf(usrdata, sortby).psm_GeneList.items(),
                                            sortdf(self.presplit, sortby).psm_GeneList.items()):
            calc = ','.split(calc)
            result = ','.split(result)
            with self.subTest(calc=calc):
                self.assertListEqual(calc, result,
                                 msg=('Error in gene list input index {}, '
                                      'result index {}').format(ix1, ix2)+ ' Sequences:' + \
                                 ' | '.join(get_sequences(ix1, usrdata,
                                                       ix2, self.presplit)))

        # test GeneCount
        for (ix1,calc), (ix2,result) in zip(sortdf(usrdata, sortby).psm_GeneCount.items(),
                                            sortdf(self.presplit, sortby).psm_GeneCount.items()):
            with self.subTest(calc=calc):
                self.assertEqual(calc, result,
                                 msg=('Error in gene count at input index {}, '
                                      'result index {}').format(ix1, ix2))
        # test ProteinList
        for (ix1,calc), (ix2,result) in zip(sortdf(usrdata, sortby).psm_ProteinList.items(),
                                            sortdf(self.presplit, sortby).psm_ProteinList.items()):
            calc = ','.split(calc)
            result = ','.split(result)
            with self.subTest(calc=calc):
                self.assertListEqual(calc, result,
                                 msg=('Error in protein list at input index {}, '
                                      'result index {}').format(ix1, ix2))
        # test ProteinCount
        for (ix1,calc), (ix2,result) in zip(sortdf(usrdata, sortby).psm_ProteinCount.items(),
                                            sortdf(self.presplit, sortby).psm_ProteinCount.items()):
            with self.subTest(calc=calc):
                self.assertEqual(calc, result,
                                 msg=('Error in protein count at input index {}, '
                                      'result index {}').format(ix1, ix2))

    def sequence_lowerAssert(self):
        """Test if sequence_lower column is equal to Sequence.lower() column"""
        usrdata = pd.DataFrame(self.presplit, copy=True)
        usrdata['sequence_lower'] = usrdata.apply(lambda x: x['Sequence'].lower(),
                                                  axis=1)
        for row in usrdata[['Sequence', 'sequence_lower']].itertuples():
            seq = row.Sequence
            with self.subTest(seq=seq):
                self.assertSequenceEqual(row.Sequence.lower(), row.sequence_lower)

    def redundant_peakAssert(self):
        """Test if redundant peaks are identified as they should be"""
        usrdata = pd.DataFrame(self.presplit, copy=True)
        usrdata.drop(labels='psm_Peak_UseFLAG', inplace=True, axis=1)
        #usrdata.to_csv('usrtest_part.tab',sep='\t')
        #self.presplit.to_csv('presplit_part.tab',sep='\t')
        #print('testing redundant peaks')
        #self.usrdata, self.presplit = update_df(self.usrdata, self.presplit, reassignments=reassign)
        #print('just updated input data')

        #self.usrdata.to_csv('usrtest_part.tab', sep='\t')
        #self.presplit.to_csv('presplit_part.tab', sep='\t')
        usrdata = pygrouper.redundant_peaks(usrdata)
        #usrdata.to_csv('usrtest1_part.tab',sep='\t')
        #tempf = open('temp.txt','w')
        for (ix1,calc), (ix2,result) in zip(usrdata.psm_Peak_UseFLAG.items(),
                                            self.presplit.psm_Peak_UseFLAG.items()):
#        for (ix1,calc), (ix2,result) in zip(usrdata.iterrows(),
#                                            self.presplit.iterrows()):
#            calc_value = calc['psm_Peak_UseFLAG']
#            result_value = result['psm_Peak_UseFLAG']
            #tempf.write('\n'+'*'*70+'\n')
            #tempf.write('*'*70+'\n')
            #tempf.write('input index : {} \n'.format(ix1))
            #tempf.write(usrdata.loc[ix1].to_string())
            #tempf.write('\n'+'*'*70+'\n')
            #tempf.write('result index : {} \n'.format(ix2))
            #tempf.write(self.presplit.loc[ix2].to_string())
            #tempf.write('\n'+'*'*70+'\n')
            #tempf.write('*'*70+'\n')
            with self.subTest(calc=calc):
                self.assertEqual(calc, result,
                                 msg=('Error in Peak UseFlag at input index {}, '
                                      'result index {}').format(ix1, ix2))

        #reassign = ['psm_Peak_UseFLAG']
        #self.usrdata = update_df(self.usrdata, self.presplit, reassignments=reassign)


    def sum_areaAssert(self):
        """Test if similar peaks are summed appropriately"""
        self.presplit.fillna(0, inplace=True)
        usrdata = pd.DataFrame(self.presplit, copy=True).fillna(0)
        usrdata.drop(labels='psm_SequenceArea', inplace=True, axis=1)
        usrdata = pygrouper.sum_area(usrdata, self.area_col)
        #usrdata.to_csv('usrtest_part.tab', sep='\t')
        #self.presplit.to_csv('presplit_part.tab', sep='\t')

        for (ix1,calc), (ix2,result) in zip(usrdata.psm_SequenceArea.items(),
                                            self.presplit.psm_SequenceArea.items()):
            with self.subTest(calc=calc):
                self.assertEqual(calc, result,
                                 msg=('Error in Sequence Area at input index {}, '
                                      'result index {}').format(ix1, ix2))

    def auc_reflagAssert(self):
        """Test if redundant peaks are properly marked for reflag"""

        self.presplit.fillna(0, inplace=True)
        usrdata = pd.DataFrame(self.presplit, copy=True).fillna(0)
        usrdata.drop(labels='AUC_reflagger', inplace=True, axis=1)
        usrdata, area_col = pygrouper.auc_reflagger(usrdata, self.area_col)
        self.assertFalse(area_col==self.area_col) # should be reassigned
        for (ix1,calc), (ix2,result) in zip(usrdata.AUC_reflagger.items(),
                                            self.presplit.AUC_reflagger.items()):
            with self.subTest(calc=calc):
                self.assertEqual(calc, result,
                                 msg=('Error in AUC Reflagger at input index {}, '
                                      'result index {}').format(ix1, ix2))

    def split_on_geneid(self):
        self.presplit.fillna(0, inplace=True)
        usrdata = pd.DataFrame(self.presplit, copy=True).fillna(0)
        usrdata.drop(labels='psm_GeneID', inplace=True, axis=1)

        self.postsplit.fillna(0, inplace=True)
        usrdata = pygrouper.split_on_geneid(usrdata)

        sortby=['Sequence', 'FirstScan', 'IonScore', 'PEP', 'psm_GeneID']
        sortdf(usrdata, sortby, inplace=True)
        sortdf(self.postsplit, sortby, inplace=True)

        #print(usrdata.psm_GeneID.values)
        #print(self.postsplit.psm_GeneID.values)

        self.assertEqual(len(usrdata), len(self.postsplit))
        #self.assertCountEqual(usrdata, self.postsplit)
        #for (ix1,calc), (ix2,result) in zip(usrdata.iterrows(), self.postsplit.iterrows()):
        #    calc_list = calc.values.tolist()
        #    result_list = result.values.tolist()
        #    with self.subTest(calc_list=calc_list):
        #        self.assertListEqual(calc_list, result_list)



class GrouperTest(GrouperTestCase):

    def setUp(self):
        """load data, runs once for each test"""
        dtype_dict = {'psm_GeneList': 'object', 'psm_ProteinList': 'object',
                      'psm_GeneID': 'object'}
        self.usrdata   = pd.read_table(os.path.join(self.path, self.inputf))
        self.refseq    = pd.read_table(os.path.join(self.path, self.refseqf))
        self.presplit  = pd.read_table(os.path.join(self.path, self.presplitf),
                                       dtype=dtype_dict) # pre-split results
        self.postsplit = pd.read_table(os.path.join(self.path, self.postsplitf),
                                       dtype=dtype_dict) # post-split results

        self.e2gdata   = pd.read_table(os.path.join(self.path, self.e2gfile))
        self.area_col  = 'PrecursorArea'
#        print('calling setup')
    #def tearDown(self)
    #""" clean up files written to disk """
    #    # tear down

    def test_seqmodi(self):
        self.seq_modiAssert()

    def test_genemap(self):
        self.gene_mapperAssert()

    def test_seqlower(self):
        self.sequence_lowerAssert()

    def test_redundant_peaks(self):
        self.redundant_peakAssert()

    #@unittest.skipIf(self.currentTest == 'simple')
    def test_sequence_area(self):
        self.sum_areaAssert()

    def test_auc_reflag(self):
        self.auc_reflagAssert()

    def test_split_on_geneid(self):
        self.split_on_geneid()



def suite():
    suite = unittest.TestSuite()

    test_parameters = get_sample_data()
    for test_parameter in test_parameters:
        testdata = test_parameters[test_parameter]
        loadedtests = unittest.TestLoader().loadTestsFromTestCase(GrouperTest)
        for t in loadedtests:
            t.currentTest = testdata.directory
            t.path        = test_parameter
            t.inputf      = testdata.inputdata
            t.refseqf     = testdata.refseq
            t.presplitf   = testdata.presplit
            t.postsplitf  = testdata.postsplit
            t.e2gfile     = testdata.genedata
        suite.addTests(loadedtests)

    return suite

def run_tests(verbose=2, hide_stdout=True, tofile=False, testrunner=None):
    """Get all test data, put it together, and run
    """
    if testrunner is None:
        testrunner = 'unittest'
    elif isinstance(testrunner, str):
        testrunner = testrunner.lower()
    if hide_stdout:
        f = open(os.devnull, 'w')
        f = open('temp_stdout.txt','w')
        sys.stdout = f # suppress all print statements from pygrouper (as well as anywhere else)
    test_suite = suite()

    testskips = [('simple', 'test_sequence_area', 'not in output yet')]
    for ts in test_suite:
        for skip in testskips:
            if skip[0] in ts.__str__() and skip[1] in ts.__str__():
                setattr(ts, 'setUp', lambda: ts.skipTest(skip[2]))
                break
    if testrunner == 'unittest':
        if tofile:
            log_file = 'test_log.txt'
            f = open(log_file, 'w')
            runner = unittest.TextTestRunner(f, verbosity=verbose)
        else:
            runner = unittest.TextTestRunner(verbosity=verbose)
        runner.run(test_suite)
    elif testrunner == 'html':
        import HTMLTestRunner
        fp = open('my_testing_report.html', 'wb')
        runner = HTMLTestRunner.HTMLTestRunner(
            stream=fp,
            title='My pygrouper unit test',
            description='Testing pygrouper'
        )
        runner.run(test_suite)
    elif testrunner == 'nose':
        import nose
        nose.run(suite=nose.suite.LazySuite(test_suite))

if __name__ == '__main__':
    run_tests(tofile=False, hide_stdout=False)
    run_tests(tofile=True)
    #result = unittest.TestResult()
    #test_suite.run(result)
