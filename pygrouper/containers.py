"""Container for each experiment, has a dataframe and metadata"""
import os
from datetime import datetime
import pandas as pd

class UserData:

    def __init__(self, datafile=None, runno=1, searchno=1, no_taxa_redistrib=0, addedby='',
                 indir = '.', outdir='.', rawfiledir='.', usedb=False, labeltype='none',
                 searchdb=None):
        self.recno = None
        self.runno = runno
        self.searchno = searchno
        self.taxonid = None
        self.quant_source = None
        self.added_by = addedby
        self.techrepno = 1 # depreciated
        self.labeltype = labeltype
        self.no_taxa_redistrib = no_taxa_redistrib
        self.usedb = usedb
        self.filtervalues = dict()
        self.indir = indir
        self.outdir = outdir
        self.rawfiledir = rawfiledir
        self.searchdb = searchdb # file name for refseq
        self.datafile = datafile
        self.df = pd.DataFrame()
        self.pipeline = None
        self.original_columns = None


    def __repr__(self):
        return '{}_{}_{}'.format(self.recno, self.runno, self.searchno)

    def __bool__(self):
        if self.datafile is not None and self.recno is not None:
            return True
        return False

    def full_path(self, in_or_out='in'):
        """returns data file with given path"""
        if in_or_out == 'in':
            mydir = self.indir
        elif in_or_out == 'out':
            mydir = self.outdir
        else:
            mydir = '.'
        return os.path.join(mydir, self.datafile or '')

    def read_csv(self, *args, **kwargs):
        """Uses pandas read_csv function to read an input file
        args and kwargs are passed to this function"""
        self.df = pd.read_csv(self.full_path(), *args, **kwargs)
        self.original_columns = self.df.columns.values
        return self

    def output_name(self, *suffix, ext='tab'):
        """generate an appropriate output file name
        returns rec_run_search_labeltype_filetype.tab"""
        suffix = '_'.join([str(ix) for ix in suffix])
        return '{!r}_{}{}.{}'.format(self,
                                     self.labeltype,
                                     '_' + suffix if suffix else '',
                                     ext
        )

    def populate_base_data(self):
        """Populate dataframe with base data prior to grouping"""
        self.df['psm_EXPRecNo'] = self.recno
        self.df['psm_EXPRunNo'] = self.runno
        self.df['psm_EXPSearchNo'] = self.searchno
        self.df['psm_EXPTechRepNo'] = self.techrepno
        self.df['psm_CreationTS'] = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
        self.df['psm_AddedBy'] = self.added_by
        # self.df['psm_TaxonID'] = self.taxonid
        #self.df['psm_GeneList'] = ''
        #self.df['psm_ProteinList'] = ''
        #self.df['psm_GeneCount'] = 0
        #self.df['psm_ProteinCount'] = 0
        #self.df['psm_HomologeneID'] = ''
        #self.df['psm_ProteinCapacity'] = ''
        # self.df['metadatainfo'] = [tuple()] * len(self.df)
        self.df['metadatainfo'] = ''
        if not 'ion_score_bins' in self.filtervalues:
            self.filtervalues['ion_score_bins'] = (10, 20, 30)
        return self

    @property
    def filterstamp(self):
        return 'is{ion_score}_qv{qvalue}_pep{pep}_idg{idg}_z{zmin}to{zmax}_mo{modi}_is_bins{ion_score_bins}'.format(**self.filtervalues)
