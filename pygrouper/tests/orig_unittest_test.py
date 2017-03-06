import sys
import unittest
import pandas as pd
from pandas.util.testing import assert_frame_equal
sys.path.append('..')
import pygrouper
pd.set_option('max_columns',99)
def get_test1(position=None):

    #df =  pd.read_table('C:/Users/saltzman/Desktop/temp/subdf2_after.txt')
    #df['sequence_lower'] = df.apply(lambda x: x['Sequence'].lower(),
    #                                axis=1)
    #df['psm_Peak_UseFLAG'] = [1, 0, 1, 1, 1, 1, 1, 0]
    #return df
    d = {'PrecursorArea':  [124010487227.5, 124010487227.5,
                             1415145952, 284121588, 259877834,
                             172547188, 0, 0],
         'PSMAmbiguity': 'unambiguous',
         'IonScore': [73.89, 62.93, 49.22, 47.89, 66.64, 73.73,
                      94.07, 89.16],
         'SpectrumFile': ['30381_1_QEP_ML021_75min_020_RP21.raw',
                          '30381_1_QEP_ML021_75min_020_RP21.raw',
                          '30381_1_QEP_ML021_75min_020_RP21.raw',
                          '30381_1_QEP_ML021_75min_020_RP21.raw',
                          '30381_1_QEP_ML021_75min_020_RP0625.raw',
                          '30381_1_QEP_ML021_75min_020_RP0930.raw',
                          '30381_1_QEP_ML021_75min_020_RP21.raw',
                          '30381_1_QEP_ML021_75min_020_RP21.raw',],
         'Charge': 3,
         'MisedCleavages': 0,
         'q_value' : 0,
         'psm_SequenceModiCount': 0,
         'psm_LabelFLAG': 0,
         'Sequence': 'VIISAPSADAPMFVMGVNHEK',
         'psm_SequenceModi': 'VIISAPSADAPMFVMGVNHEK',
         'FirstScan': [18060, 18235, 18555, 18715, 25605,
                       25717, 18390, 18890],
         'RTmin': [44.744115, 45.130631, 45.879951, 46.255509,
                   59.060033, 59.522358, 45.507751, 46.632535],
         'PEP': [3.84*10**-10, 6.17*10**-8, 1.62*10**-7,
                 5.90*10**-10, 1.50*10**-12, 1.16*10**-7,
                 7.05*10**-9, 3.16*10**-11],}


    df = pd.DataFrame(d)
    if position is None:
        return df
    df['psm_GeneList'] = '14433,100042025,2597'
    df['psm_GeneCount'] = 3
    # not currently tested, but have never seen an issue
    df['psm_PSM_IDG'] = df.apply(lambda x:
                                 pygrouper.IDG_picker(x['IonScore'],
                                            x['q_value']), axis=1) 
    
    if position == 'after_matched':
        return df
    df['sequence_lower'] = df.apply(lambda x: x['Sequence'].lower(),
                                    axis=1)
    if position == 'after_seq_lower':
        return df

    df['psm_Peak_UseFLAG'] = [1, 0, 1, 1, 1, 1, 1, 0]
    if position == 'after_redundant':
        return df

    df['psm_SequenceArea'] = 126142179789.5
    if position == 'after_areasum':
        return df
    df['AUC_reflagger'] = [1, 0, 0, 0, 0, 0, 0, 0]
    if position == 'after_auc_reflag':
        return df

    gid1 = pd.DataFrame(df, copy=True)
    gid2 = pd.DataFrame(df, copy=True)
    gid3 = pd.DataFrame(df, copy=True)
    gid1['psm_GeneID'] = '2597'
    gid2['psm_GeneID'] = '14433'
    gid3['psm_GeneID'] = '100042025'
    df = pd.concat([gid1, gid2, gid3], ignore_index=True)
    #sort_by = ['psm_GeneID', 'FirstScan', 'AUC_reflagger', 'psm_Peak_UseFLAG']
    #print(df.sort_values(sort_by))
    if position == 'after_gid_split':
        return df
    df['psm_AUC_useflag'] = [1, 0, 0, 0, 0, 0, 0, 0, 1,
                             0, 0, 0, 0, 0, 0, 0, 1, 0,
                             0, 0, 0, 0, 0, 0, ]
    df['psm_PSM_useflag'] = 1
    if position == 'after_flagging':
        return df
    
    print('*'*10,'','You got too far!','','*'*10, sep='\n')

def get_tempdf(usrdata, label=1):
    temp_df = usrdata[(usrdata['psm_LabelFLAG'] == label) &
                      (usrdata['psm_AUC_useflag'] == 1) &
                      (usrdata['psm_GeneCount'] > 0)]  # should keep WL's
    #peptides too (if they have a mass)
    return temp_df

class GrouperTestCase(unittest.TestCase):


    def sequence_lower(self, usrdata):
        """Test if sequence_lower column is equal to Sequence.lower() column"""
        usrdata['sequence_lower'] = usrdata.apply(lambda x: x['Sequence'].lower(),
                                                  axis=1)
        for row in usrdata[['Sequence', 'sequence_lower']].itertuples():
            seq = row.Sequence
            with self.subTest(seq=seq):
                self.assertSequenceEqual(row.Sequence.lower(), row.sequence_lower)

    def redundant_peaks(self, usrdata, result):
        """Test if redundant peaks are identified as they should be"""
        usrdata = pygrouper.redundant_peaks(usrdata)
        for calc, actual in zip(usrdata['psm_Peak_UseFLAG'].values,
                                result):
            with self.subTest(calc=calc):
                self.assertEqual(calc, actual)

    def sum_areaAssert(self, usrdata, area_col, result):
        """Test if similar peaks are summed appropriately"""
        usrdata = pygrouper.sum_area(usrdata, area_col)
        for calc, actual in zip(usrdata['psm_SequenceArea'].values,
                                result):
            with self.subTest(calc=calc):
                self.assertEqual(calc, actual)

    def auc_reflag(self, usrdata, area_col, result):
        """Test if redundant peaks are properly marked for reflag"""
        usrdata, area_col = pygrouper.auc_reflagger(usrdata, area_col)
        for calc, actual in zip(usrdata['AUC_reflagger'].values,
                                result):
            with self.subTest(calc=calc):
                self.assertEqual(calc, actual)

    def geneid_split(self, usrdata, result):
        """Test for proper splitting on geneid
        Note: GeneID column should be an object (text)
        Note: make sure you sort and drop indexes of both before comparing"""
        usrdata = pygrouper.split_on_geneid(usrdata)
        sort_by = ['psm_GeneID', 'FirstScan', 'AUC_reflagger', 'psm_Peak_UseFLAG']
        self.assertTrue(usrdata.sort_values(sort_by).reset_index(drop=True).\
                        equals(result.sort_values(sort_by).reset_index(drop=True)))

    def correct_flagAssert(self, usrdata, result):
        """Test for proper AUC and PSM flag assignments"""
        FilterValues = {'Filter_IS': 7, 'Filter_qV': 0.05, 'Filter_PEP': 'all',
                        'Filter_IDG': 'all', 'Filter_Z_min': 2, 'Filter_Z_max': 4,
                        'Filter_Modi': 3}  # defaults
        sort_by = ['psm_GeneID', 'FirstScan', 'AUC_reflagger', 'psm_Peak_UseFLAG']
        new =['psm_AUC_useflag', 'psm_PSM_useflag']
        usrdata['psm_AUC_useflag'], usrdata['psm_PSM_useflag'] = \
        list(zip(*usrdata.apply(pygrouper.AUC_PSM_flagger, args=(FilterValues,), axis=1)))
        for (calcs, actuals) in zip(usrdata.sort_values(sort_by).reset_index(drop=True)[new].values,
                                    result.sort_values(sort_by).reset_index(drop=True)[new].values):
            with self.subTest(calcs=calcs):
                self.assertListEqual(calcs.tolist(), actuals.tolist())

    def gene_sets(self, usrdata, label=1):
        temp_df = get_tempdf(usrdata,label=label)
        genes_df = pd.DataFrame({'gene_GeneID':
                                 list(set(temp_df['psm_GeneID']))})

        self.assertSetEqual(set(genes_df.gene_GeneID.tolist()),
                            set(temp_df.psm_GeneID.tolist()))

        def peptide_to_genes(self, usrdata):
            pass

        #print(usrdata.head())
        #print()
        #print(result.head())
#        self.assertTrue(assert_frame_equal(usrdata.sort_values(sort_by).reset_index(drop=True),
#                                           result.sort_values(sort_by).reset_index(drop=True)))
                                           
        #self.assertTrue(usrdata.sort_values(sort_by).reset_index(drop=True).\
        #                equals(result.sort_values(sort_by).reset_index(drop=True)))


class GrouperTest(GrouperTestCase):
    #def setUp(self):
    #    pass
    #def test_easy(self):
    #    pass
    #def test_hard(self):
    #    pass
    def test1(self):
        #test1
        data = get_test1(position='after_matched')
        self.sequence_lower(data)

        #test2
        data = get_test1(position='after_seq_lower')
        data_results = get_test1(position='after_redundant')
        result = data_results['psm_Peak_UseFLAG'].values
        self.redundant_peaks(data, result)
        del data, data_results
        #test3
#'after_redundant'
        data = get_test1(position='after_redundant')
        data_results = get_test1(position='after_areasum')
        result = data_results['psm_SequenceArea'].values
        self.sum_areaAssert(data, 'PrecursorArea', result)

        data = data_results
        data_results = get_test1('after_auc_reflag')
        result = data_results['AUC_reflagger'].values
        self.auc_reflag(data, 'PrecursorArea', result)

        data = data_results
        data_results = get_test1('after_gid_split')
        self.geneid_split(data, data_results)

        data = data_results
        data_results = get_test1('after_flagging')
        self.correct_flagAssert(data, data_results)

        data = data_results
        self.gene_sets(data)

if __name__ == '__main__':
    unittest.main()
