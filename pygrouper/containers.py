"""Container for each experiment, has a dataframe and metadata"""
import os
from datetime import datetime
import traceback

import pandas as pd


from . import _version

class UserData:

    def __init__(self, recno=None, datafile=None, runno=1, searchno=1, no_taxa_redistrib=0,
                 addedby='', indir = '.', outdir='.', rawfiledir='.',
                 labeltype='none', quant_source=None, phospho=False,
                 searchdb=None, taxonid=None, miscuts=2):
        if recno is None:
            raise ValueError('Must supply record number (recno)')
        self.recno = recno
        self.runno = runno
        self.searchno = searchno
        self.taxonid = taxonid
        self.added_by = addedby
        self.labeltype = labeltype
        self.no_taxa_redistrib = no_taxa_redistrib
        self.filtervalues = dict()
        self.indir = indir
        self.outdir = outdir
        self.rawfiledir = rawfiledir
        self.searchdb = searchdb # file name for refseq
        self.datafile = datafile
        self.df = pd.DataFrame()
        self.pipeline = None
        self.original_columns = None
        self.LOGFILE = os.path.join(outdir, self.output_name(ext='log'))
        self._LOGSTACK = list()
        self.EXIT_CODE = 0
        self.ERROR = None
        self.taxon_ratio_totals = dict()
        self.miscuts = miscuts
        self.phospho = phospho

        with open(self.LOGFILE, 'w') as f:
            f.write('{} PyGrouper {}'.format(datetime.now(), _version.__version__))

    @property
    def taxon_miscut_id(self):
        return hash(self.taxonid) + hash(self.miscuts)


    def __repr__(self):
        return '{}_{}_{}'.format(self.recno, self.runno, self.searchno)

    def __bool__(self):
        if self.datafile is not None and self.recno is not None:
            return True
        return False

    def to_log(self, message):
        if self._LOGSTACK:  # flush
            messages = self._LOGSTACK + (messages,)
        else:
            messages = (message,)
        with open(self.LOGFILE, 'w+') as f:
            for message in messages:
                f.write(message)
                # f.write(sep)
                f.write('\n')

    def to_logq(self, message):
        self._LOGSTACK.append(message+'\n')
        return self

    def flush_log(self):
        if self._LOGSTACK:
            stack, self._LOGSTACK = self._LOGSTACK, list()
            self.to_log('\n'.join(stack))
        return self


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
        try:
            self.df = pd.read_csv(self.full_path(), *args, **kwargs)
            self.original_columns = self.df.columns.values
        except Exception as e:
            # self.to_log(''.join(traceback.format_exc()))
            self.to_log(traceback.format_exc())
            self.ERROR = traceback.format_exc()
            self.EXIT_CODE = 1
            return 1
        if len(self.df) == 0:
            self.EXIT_CODE = 1
            return 2
        return 0

    def output_name(self, suffix=None, ext='tab'):
        """generate an appropriate output file name
        returns rec_run_search_labeltype_filetype.tab"""
        # suffix = '_'.join([str(ix) for ix in suffix])
        return '{!r}_{}{}.{}'.format(self,
                                     self.labeltype,
                                     '_' + suffix if suffix else '',
                                     ext
        )

    def populate_base_data(self):
        """Populate dataframe with base data prior to grouping"""

        self.categorical_assign('EXPRecNo', self.recno)
        self.categorical_assign('EXPRunNo', self.runno)
        self.categorical_assign('EXPSearchNo', self.searchno)
        self.categorical_assign('CreationTS', datetime.now().strftime("%m/%d/%Y) %H:%M:%S"))
        self.categorical_assign('AddedBy', self.added_by)
        # self.categorical_assign('metadatainfo', '')  # not sure if this is okay

        # self.df['EXPRecNo'] = self._categorical_assign(self.recno)
        # self.df['EXPRunNo'] = self._categorical_assign(self.runno)
        # self.df['EXPSearchNo'] = self._categorical_assign(self.searchno)
        # self.df['CreationTS'] = self._categorical_assign(datetime.now().strftime("%m/%d/%Y) %H:%M:%S"))
        # self.df['AddedBy'] = self._categorical_assign(self.added_by)

        # self.df['psm_EXPTechRepNo'] = self.techrepno
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
        s = 'is{ion_score}_qv{qvalue}_pep{pep}_idg{idg}_z{zmin}to{zmax}_mo{modi}_is_bins{ion_score_bins}'.format(**self.filtervalues)
        if self.phospho:
            s += '_phospho_only'
        return s

    def categorical_assign(self, name, value, **kwargs):
        """
        Assign a static value to a new column.
        Saves memory by using pandas Categorical dtype.
        :kwargs: passed to pd.Series.astype
        """
        self.df[name] = value
        self.df[name] = self.df[name].astype('category', **kwargs)
        return self
