#===============================================================================#
# PyGrouper - Alex Saltzman
from __future__ import print_function

import re, os, sys, time
import itertools
import json
import logging
from time import sleep
from collections import defaultdict
from functools import partial
from math import ceil
from warnings import warn
import six
if six.PY3:
    from configparser import ConfigParser
elif six.PY2:
    from ConfigParser import ConfigParser
from itertools import repeat
import traceback
import multiprocessing

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

from RefProtDB.utils import fasta_dict_from_file

from . import _version
from .subfuncts import *

# from ._orig_code import timed
pd.set_option(
    "display.width", 170,
    "display.max_columns", 500,
)
__author__ = 'Alexander B. Saltzman'
__copyright__ = _version.__copyright__
__credits__ = ['Alexander B. Saltzman', 'Anna Malovannaya']
__license__ = 'MIT'
__version__ = _version.__version__
__maintainer__ = 'Alexander B. Saltzman'
__email__ = 'saltzman@bcm.edu'
program_title = 'Pygrouper v{}'.format(__version__)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
logfilename = program_title.replace(' ', '_') + '.log'
logging.basicConfig(filename=logfilename, level=logging.DEBUG)
logging.info('{}: Initiating {}'.format(datetime.now(), program_title))

SEP = ';'

labelflag = {'none': 0,  # hard coded number IDs for labels
             'TMT_126': 1260,
             'TMT_127_C': 1270,
             'TMT_127_N': 1271,
             'TMT_128_C': 1280,
             'TMT_128_N': 1281,
             'TMT_129_C': 1290,
             'TMT_129_N': 1291,
             'TMT_130_C': 1300,
             'TMT_130_N': 1301,
             'TMT_131': 1310,
             'iTRAQ_114': 113,
             'iTRAQ_114': 114,
             'iTRAQ_115': 115,
             'iTRAQ_116': 116,
             'iTRAQ_117': 117,
             'iTRAQ_118': 118,
             'iTRAQ_119': 119,
             'iTRAQ_121': 121,
}

E2G_COLS = ['EXPRecNo', 'EXPRunNo', 'EXPSearchNo',
            'EXPLabelFLAG', 'AddedBy',
            'CreationTS', 'ModificationTS', 'TaxonID',
            'GeneID', 'ProteinGIs', 'ProteinRefs',
            'HIDs', 'Description',
            'IDSet', 'IDGroup', 'IDGroup_u2g',
            'SRA',
            'GPGroup', 'GPGroups_All', 'PSMs',
            'PSMs_u2g', 'PeptidePrint', 'PeptideCount',
            'PeptideCount_u2g', 'PeptideCount_S',
            'PeptideCount_S_u2g',
            'AreaSum_u2g_0', 'AreaSum_u2g_all',
            'AreaSum_max', 'AreaSum_dstrAdj',
            'GeneCapacity', 'iBAQ_dstrAdj']

DATA_COLS = ['EXPRecNo', 'EXPRunNo', 'EXPSearchNo',
             'Sequence', 'PSMAmbiguity', 'Modifications',
             'ActivationType', 'DeltaScore', 'DeltaCn',
             'Rank', 'SearchEngineRank', 'PrecursorArea',
             'q_value', 'PEP', 'IonScore',
             'MissedCleavages', 'IsolationInterference', 'IonInjectTime',
             'Charge', 'mzDa', 'MHDa',
             'DeltaMassDa', 'DeltaMassPPM', 'RTmin',
             'FirstScan', 'LastScan', 'MSOrder', 'MatchedIons',
             'SpectrumFile', 'AddedBy',
             'oriFLAG',
             'CreationTS', 'ModificationTS', 'GeneID',
             'GeneIDs_All', 'GeneIDCount_All',
             'ProteinGIs',
             'ProteinGIs_All', 'ProteinGICount_All',
             'ProteinRefs',
             'ProteinRefs_All', 'ProteinRefCount_All',
             'HIDs', 'HIDCount_All',
             'TaxonID', 'TaxonIDs_All', 'TaxonIDCount_All',
             'PSM_IDG', 'SequenceModi',
             'SequenceModiCount', 'LabelFLAG',
             'PeptRank', 'AUC_UseFLAG', 'PSM_UseFLAG',
             'Peak_UseFLAG', 'SequenceArea', 'PrecursorArea_split',
             # 'RazorArea',
             'PrecursorArea_dstrAdj']
_EXTRA_COLS = ['LastScan', 'MSOrder', 'MatchedIons']  # these columns are not required to be in the output data columns

try:
    from PIL import Image, ImageFont, ImageDraw
    imagetitle = True
except ImportError:
    imagetitle = False


def _apply_df(input_args):
    df, func, i, func_args, kwargs = input_args
    return i, df.apply(func, args=(func_args), **kwargs)

def apply_by_multiprocessing(df, func, workers=1, func_args=None, **kwargs):
    """
    Spawns multiple processes if has os.fork and workers > 1
    """
    if func_args is None:
        func_args = tuple()

    if workers == 1 or not hasattr(os, 'fork'):
        result = _apply_df((df, func, 0, func_args, kwargs,))
        return result[1]


    pool = multiprocessing.Pool(processes=workers)
    with multiprocessing.Pool(processes=workers) as pool:

        result = pool.map(_apply_df, [(d, func, i, func_args, kwargs,)
                                      for i, d in enumerate(np.array_split(df, workers))]
        )
        # pool.close()

    result = sorted(result, key=lambda x: x[0])

    return pd.concat([x[1] for x in result])


def quick_save(df,name='df_snapshot.p', path=None, q=False):
    import pickle
    #import RefSeqInfo

    if path:
        name = path+name
    #df.to_csv('test_matched.tab', index=False, sep='\t')
    pickle.dump(df, open(name, 'wb'))
    print('Pickling...')
    if q:
        print('Exiting prematurely')
        sys.exit(0)

def _get_rawfile_info(path, spectraf):
    if path is None:
        path = '.'
    if not os.path.isdir(path):
        return ('not found, check rawfile path', 'not found')
    for f in os.listdir(path):
        if f == spectraf:
            rawfile = os.path.abspath(os.path.join(path,f))
            break
    else:
        return ('not found', 'not found')

    fstats   = os.stat(rawfile)
    mod_date = datetime.fromtimestamp(fstats.st_mtime).strftime("%m/%d/%Y %H:%M:%S")
    size     = byte_formatter(fstats.st_size)
    return (size, mod_date)

def _spectra_summary(spectraf, data):
    """ Calculates metadata per spectra file.
    The return order is as follows:

    -minimum RT_min
    -maximum RT_min
    -min IonScore
    -max IonScore
    -min q_value
    -max q_value
    -min PEP
    -max PEP
    -min Area (precursor, exculding zeros)
    -max Area
    -PSM Count
    -median DeltaMassPPM
    """
    data      = data[data.SpectrumFile==spectraf]

    RT_min       = data.RTmin.min()
    RT_max       = data.RTmin.max()
    IonScore_min = data.IonScore.min()
    IonScore_max = data.IonScore.max()
    q_min        = data.q_value.min()
    q_max        = data.q_value.max()
    PEP_min      = data.PEP.min()
    PEP_max      = data.PEP.max()
    area_min     = data[data.PrecursorArea!=0].PrecursorArea.min()
    area_max     = data.PrecursorArea.max()
    count    = len(data[data.PSM_UseFLAG==1])
    dmass_median = data.DeltaMassPPM.median()
    return(RT_min, RT_max, IonScore_min, IonScore_max, q_min, q_max,
           PEP_min, PEP_max, area_min, area_max, count, dmass_median)

def spectra_summary(usrdata):
    """Summaries the spectral files included in an analysis.

    Args:
        usrdata: a UserData instance with the data loaded


    Returns:
        A pandas DataFrame with relevant columns, ready for export

        if the raw files cannot be found at usrdata.rawfiledir,
        then 'not found' is returned for those columns
    """
    msfdata = pd.DataFrame()
    # msfdata['RawFileName']    = list(set(usrdata.df.SpectrumFile.tolist()))
    msfdata['RawFileName']    = sorted(usrdata.df.SpectrumFile.unique())
    msfdata['EXPRecNo']       = usrdata.recno
    msfdata['EXPRunNo']       = usrdata.runno
    msfdata['EXPSearchNo']    = usrdata.searchno
    msfdata['AddedBy']        = usrdata.added_by
    msfdata['CreationTS']     = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
    msfdata['ModificationTS'] = datetime.now().strftime("%m/%d/%Y %H:%M:%S")

    summary_info = msfdata.apply(lambda x:
                        _spectra_summary(x['RawFileName'],
                                        usrdata.df),
                        axis=1)

    (msfdata['RTmin_min'], msfdata['RTmin_max'], msfdata['IonScore_min'],
     msfdata['IonScore_max'], msfdata['qValue_min'], msfdata['qValue_max'],
     msfdata['PEP_min'], msfdata['PEP_max'], msfdata['Area_min'],
     msfdata['Area_max'], msfdata['PSMCount'],
     msfdata['DeltaMassPPM_med']) = zip(*summary_info)

    rawfile_info = msfdata.apply(lambda x:
                                    _get_rawfile_info(usrdata.rawfiledir,
                                                      x['RawFileName']),
                                 axis=1)

    msfdata['RawFileSize'], msfdata['RawFileTS'] = zip(*rawfile_info)

    return msfdata

def get_gid_ignore_list(inputfile):
    """Input a file with a list of geneids to ignore when normalizing across taxa
    Each line should have 1 geneid on it.
    Use '#' at the start of the line for comments
    Output a list of geneids to ignore.
    """
    # Don't convert GIDs to ints,
    # GIDs are not ints for the input data
    return [x.strip() for x in open(inputfile, 'r') if
            not x.strip().startswith('#')]

def _format_peptideinfo(row):

    if len(row) == 0:
        return ('', 0, '', 0, '', 0, '', 0, '', 0, ())
    result = (
        # ','.join(row['GeneID'].dropna().unique()),
        SEP.join(str(x) for x in set(row['geneid'])),
        row['geneid'].replace('', np.nan).nunique(dropna=True),

        # ','.join(row['TaxonID'].dropna().unique()),
        SEP.join(str(x) for x in set(row['taxon'])),
        row['taxon'].replace('', np.nan).nunique(dropna=True),

        # ','.join(row['ProteinGI'].dropna().unique()),
        SEP.join(str(x) for x in set(row['gi'])),
        row['gi'].replace('', np.nan).nunique(dropna=True),

        SEP.join(str(x) for x in set(row['ref'])),
        row['ref'].replace('', np.nan).nunique(dropna=True),

        # ','.join(row['HomologeneID'].dropna().unique()),
        SEP.join(str(x) for x in set(row['homologene'])),
        row['homologene'].replace('', np.nan).nunique(dropna=True),

        tuple(row['capacity']),
          # row['capacity'].mean(),

    )
    return result

def _extract_peptideinfo(row, database):
    return _format_peptideinfo(database.loc[row])

def extract_peptideinfo(usrdata, database):
    filter_int = partial(filter, lambda x : x.isdigit())
    to_int = partial(map, int)
    ixs = (usrdata.df.metadatainfo.str.strip('|')
           .str.split('|')
           .apply(filter_int)
           .apply(to_int)
           .apply(list)
           # .apply(pd.Series)
           # .stack()
           # .to_frame()
           )

    # info = ixs.apply(lambda x : _format_peptideinfo(database.loc[x])).apply(pd.Series)
    info = ixs.pipe(apply_by_multiprocessing,
                    _extract_peptideinfo,
                    func_args=(database,),
                    workers=WORKERS,
    ).apply(pd.Series)

    info.columns = ['GeneIDs_All', 'GeneIDCount_All', 'TaxonIDs_All', 'TaxonIDCount_All', 'ProteinGIs_All',
                    'ProteinGICount_All', 'ProteinRefs_All', 'ProteinRefCount_All', 'HIDs', 'HIDCount_All',
                    'GeneCapacities']
    for col in ('TaxonIDs_All', 'TaxonIDCount_All', 'ProteinGIs_All', 'ProteinGICount_All',
                'ProteinRefs_All', 'ProteinRefCount_All', 'HIDs', 'HIDCount_All'):
        info[col] = info[col].astype('category')

    usrdata.df = usrdata.df.join(info)
    # (usrdata.df['GeneIDs_All'],
    #  usrdata.df['GeneIDCount_All'],
    #  usrdata.df['TaxonIDs_All'],
    #  usrdata.df['TaxonIDCount_All'],
    #  usrdata.df['ProteinGIs_All'],
    #  usrdata.df['ProteinGICount_All'],
    #  usrdata.df['ProteinRefs_All'],
    #  usrdata.df['ProteinRefCount_All'],
    #  usrdata.df['HIDs'],
    #  usrdata.df['HIDCount_All'],
    #  usrdata.df['GeneCapacities']) = zip(*info)
    # usrdata.df['TaxonIDs_All'] = usrdata.df['TaxonIDs_All'].dropna().astype(str)
    # usrdata.df['HIDs'] = usrdata.df['HIDs'].fillna('')

    return 0

def gene_mapper(df, other_col=None):

    if other_col is None or other_col not in df.columns:
        raise ValueError("Must specify other column")
    groupdf = (df[['geneid', other_col]]
               .drop_duplicates()
               .groupby('geneid')
    )

    # d = {k: SEP.join(filter(None, str(v))) for k, v in groupdf[other_col]}
    d = {k: SEP.join(filter(None, map(str, v))) for k, v in groupdf[other_col]}

    return d

def gene_taxon_mapper(df):
    """Returns a dictionary with mapping:
    gene -> taxon
    Input is the metadata extracted previously"""
    return gene_mapper(df, 'taxon')

def gene_symbol_mapper(df):
    """Returns a dictionary with mapping:
    gene -> taxon
    Input is the metadata extracted previously"""
    return gene_mapper(df, 'symbol')

def gene_desc_mapper(df):
    """Returns a dictionary with mapping:
    gene -> taxon
    Input is the metadata extracted previously"""
    return gene_mapper(df, 'description')

def gene_hid_mapper(df):
    """Returns a dictionary with mapping:
    gene -> taxon
    Input is the metadata extracted previously"""
    return gene_mapper(df, 'homologene')

def gene_protgi_mapper(df):
    """Returns a dictionary with mapping:
    gene -> taxon
    Input is the metadata extracted previously"""
    return gene_mapper(df, 'gi')

def gene_protref_mapper(df):
    """Returns a dictionary with mapping:
    gene -> taxon
    Input is the metadata extracted previously"""
    return gene_mapper(df, 'ref')




def assign_IDG(df, filtervalues=None):
    filtervalues = filtervalues or dict()
    ion_score_bins = filtervalues.get('ion_score_bins', (10, 20, 30))
    df['PSM_IDG'] = pd.cut(df['IonScore'],
                               # bins=(0, *ion_score_bins, np.inf),
                               bins=(0,) + tuple(ion_score_bins) + (np.inf,),
                               labels=[7, 5, 3, 1], include_lowest=True,
                               right=False).astype('int')
    df.loc[ df['q_value'] > .01, 'PSM_IDG' ] += 1
    return df

def make_seqlower(usrdata, col='Sequence'):
    """Make a new column called sequence_lower from a DataFrame"""
    usrdata['sequence_lower'] = usrdata[col].str.lower()
    return

def peptidome_matcher(usrdata, ref_dict):

    if not ref_dict:
        return usrdata
    ref_dict_filtered = ref_dict
    pmap = partial(map, str)
    result = (usrdata.Sequence.str.upper().map(ref_dict)
              .fillna('')
              .map(pmap)
              .map('|'.join)
              .add('|')
              )
    usrdata['metadatainfo'] += result
    return usrdata

def redundant_peaks(usrdata):
    """ Remove redundant, often ambiguous peaks by keeping the peak
    with the highest ion score"""
    peaks = usrdata.sort_values(by='IonScore', ascending=False).\
            drop_duplicates(subset=['SpectrumFile','SequenceModi', 'Charge', 'PrecursorArea'])
    peaks['Peak_UseFLAG'] = 1
    # peaks['Peak_UseFLAG'] = True
    usrdata = usrdata.join(peaks['Peak_UseFLAG'])
    usrdata['Peak_UseFLAG'] = usrdata.Peak_UseFLAG.fillna(0).astype(np.int8)
    # usrdata['Peak_UseFLAG'] = usrdata.Peak_UseFLAG.fillna(False)
    print('Redundant peak areas removed : ', len(usrdata)-len(peaks))
    return usrdata

def sum_area(df):
    """Sum the area of similar peaks
    New column SequenceArea is created
    """
    df['Sequence_set'] = df['Sequence'].apply(lambda x: tuple(set(list(x))))
    summed_area = (df.query('Peak_UseFLAG==1')
                   .filter(items=['SequenceModi', 'Charge', 'PrecursorArea'])
                   .groupby(['SequenceModi', 'Charge'])
                   .agg({'PrecursorArea': 'sum'})
                   .reset_index()
                   .rename(columns={'PrecursorArea': 'SequenceArea'})
    )
    df = df.merge(summed_area, how='left', on=['SequenceModi', 'Charge'])
    return df

def auc_reflagger(df):
    """Remove duplicate sequence areas
    """
    #usrdata['Sequence_set'] = usrdata['Sequence'].apply(lambda x: tuple(set(list(x))))
    no_dups = (df.sort_values(by=['SequenceModi', 'Charge', 'SequenceArea',
                                  'PSM_IDG', 'IonScore', 'PEP', 'q_value'],
                              ascending=[1,1,0,1,0,1,1])
               .drop_duplicates(subset=['SequenceArea', 'Charge', 'SequenceModi',])
               .assign(AUC_reflagger = True)
    )
    df = (df.join(no_dups[['AUC_reflagger']])
          .assign(AUC_reflagger = lambda x: (x['AUC_reflagger']
                                             .fillna(0)
                                             .astype(np.int8)))
    )
    return df

def export_metadata(program_title='version',usrdata=None, matched_psms=0, unmatched_psms=0,
                    usrfile='file', taxon_totals=dict(), outname=None, outpath='.', **kwargs):
    """Update iSPEC database with some metadata information
    """
    print('{} | Exporting metadata'.format(time.ctime()))
    #print('Number of matched psms : ', matched_psms)
    d = dict(
        version=program_title,
        searchdb=usrdata.searchdb,
        filterstamp=usrdata.filterstamp,
        matched_psms=matched_psms,
        unmatched_psms=unmatched_psms,
        inputname=usrdata.datafile,
        hu=taxon_totals.get('9606', 0),
        mou=taxon_totals.get('10090', 0),
        gg=taxon_totals.get('9031', 0),
        recno=usrdata.recno,
        runno=usrdata.runno,
        searchno=usrdata.searchno
    )
    with open(os.path.join(outpath, outname), 'w') as f:
        json.dump(d, f)

def split_on_geneid(df):
    """Duplicate psms based on geneids. Areas of each psm is recalculated based on
    unique peptides unique for its particular geneid later.
    """
    oriflag = lambda x: 1 if x[-1] == 0 else 0
    glstsplitter = (df['GeneIDs_All'].str.split(SEP)
                    .apply(pd.Series, 1).stack()
                    .to_frame(name='GeneID')
                    .assign(oriFLAG= lambda x: x.index.map(oriflag))
                    )

    glstsplitter.index = glstsplitter.index.droplevel(-1)  # get rid of
                                                           # multi-index
    df = (df.join(glstsplitter)
          .reset_index())
    df['GeneID'] = df.GeneID.astype(int)
    return df

def rank_peptides(df, area_col, ranks_only=False):
    """Rank peptides here
    area_col is sequence area_calculator
    ranks_only returns just the ranks column. This does not reset the original index
    """
    df = df.sort_values(by=['GeneID', area_col,
                            'SequenceModi',
                            'Charge', 'PSM_IDG', 'IonScore', 'PEP',
                            'q_value'],
                        ascending=[1, 0, 0, 1, 1, 0, 1, 1])
    if not ranks_only:  # don't reset index for just the ranks
        df.reset_index(inplace=True)  # drop=True ?
    df.Modifications.fillna('', inplace=True)  # must do this to compare nans
    df[area_col].fillna(0, inplace=True)  # must do this to compare
    #nans
    ranks = (df[ (df.AUC_UseFLAG == 1) &
                 (df.PSM_UseFLAG == 1) &
                 (df.Peak_UseFLAG == 1) ]
             .groupby(['GeneID', 'LabelFLAG'])
             .cumcount() + 1)  # add 1 to start the peptide rank at 1, not 0
    ranks.name = 'PeptRank'
    if ranks_only:
        return ranks
    df = df.join(ranks)
    return df


def flag_AUC_PSM(df, fv, contaminant_label='__CONTAMINANT__', phospho=False):

    if fv['pep'] =='all' : fv['pep'] = float('inf')
    if fv['idg'] =='all' : fv['idg'] = float('inf')
    df['AUC_UseFLAG'] = 1
    df['PSM_UseFLAG'] = 1
    df.loc[(df['Charge'] < fv['zmin']) | (df['Charge'] > fv['zmax']),
           ['AUC_UseFLAG', 'PSM_UseFLAG']] = 0

    df.loc[df['SequenceModiCount'] > fv['modi'],
           ['AUC_UseFLAG', 'PSM_UseFLAG']] = 0

    df.loc[(df['IonScore'].isnull() & df['q_value'].isnull()),
           ['AUC_UseFLAG', 'PSM_UseFLAG']] = 1, 0

    df.loc[df['IonScore'] < fv['ion_score'],
           ['AUC_UseFLAG', 'PSM_UseFLAG']] = 0

    df.loc[df['q_value'] > fv['qvalue'],
           ['AUC_UseFLAG', 'PSM_UseFLAG']] = 0

    df.loc[df['PEP'] > fv['pep'],
           ['AUC_UseFLAG', 'PSM_UseFLAG']] = 0
    df.loc[df['PSM_IDG'] > fv['idg'],
           ['AUC_UseFLAG', 'PSM_UseFLAG']] = 0

    df.loc[(df['Peak_UseFLAG'] == 0) & (df['PSMAmbiguity'].str.lower()=='unambiguous'),
           ['AUC_UseFLAG', 'PSM_UseFLAG']] = 1
    df.loc[(df['Peak_UseFLAG'] == 0) & (df['PSMAmbiguity'].str.lower()!='unambiguous'),
           ['AUC_UseFLAG', 'PSM_UseFLAG']] = 0

    df.loc[ df['AUC_reflagger'] == 0, 'AUC_UseFLAG'] = 0

    df.loc[ df['GeneIDs_All'].str.contains(contaminant_label), ['AUC_UseFLAG', 'PSM_UseFLAG'] ] = 0, 0

    if phospho:
        df.loc[ ~df['SequenceModi'].str.contains('pho', case=False), ['AUC_UseFLAG', 'PSM_UseFLAG'] ] = 0, 0

    return df

def gene_taxon_map(usrdata, gene_taxon_dict):
    """make 'gene_taxon_map' column per row which displays taxon for given gene"""
    usrdata['TaxonID'] = usrdata['GeneID'].map(gene_taxon_dict)
    return

def get_all_taxons(taxonidlist):
    """Return a set of all taxonids from
    usrdata.TaxonIDList"""
    taxon_ids = set(SEP.join(x for x in taxonidlist
                             if x.strip()).split(SEP))
    return taxon_ids

def multi_taxon_splitter(taxon_ids, usrdata, gid_ignore_list, area_col):
    """Plugin for multiple taxons
    Returns a dictionary with the totals for each detected taxon"""
    taxon_totals = dict()
    for taxon in taxon_ids:
        #all_others = [x for x in taxon_ids if x != taxon]
        uniq_taxon = usrdata[
            #(usrdata._data_tTaxonIDList.str.contains(taxon)) &
            #(~usrdata._data_tTaxonIDList.str.contains('|'.join(all_others)))&
            (usrdata['AUC_UseFLAG'] == 1) &
            (usrdata['TaxonID'] == str(taxon)) &
            (usrdata['TaxonIDCount_All'] == 1) &
            (~usrdata['GeneID'].isin(gid_ignore_list))
        ]
        taxon_totals[taxon] = (uniq_taxon[area_col] / uniq_taxon['GeneIDCount_All']).sum()

    tot_unique = sum(taxon_totals.values())  #sum of all unique
        # now compute ratio:
    for taxon in taxon_ids:
        taxon = str(taxon)
        try:
            percentage = taxon_totals[taxon] / tot_unique
        except ZeroDivisionError:
            warn("""This file has multiple taxa but no unique to taxa peptides.
            Please check this experiment
            """)
            percentage = 1
        taxon_totals[taxon] = percentage
        print(taxon, ' ratio : ', taxon_totals[taxon])
        #logfile.write('{} ratio : {}\n'.format(taxon, taxon_totals[taxon]))
    return taxon_totals

def create_df(inputdf, label, inputcol='GeneID'):
    """Create and return a DataFrame with gene/protein information from the input
    peptide DataFrame"""
    return pd.DataFrame({'GeneID':
                         list(set(inputdf[inputcol])),
                         'EXPLabelFLAG': labelflag.get(label, 0)})

def select_good_peptides(usrdata, labelix):
    """Selects peptides of a given label with the correct flag and at least one genecount
    The LabelFLAG is set here for TMT/iTRAQ/SILAC data.
    """
    temp_df = usrdata[(usrdata['PSM_UseFLAG'] == 1) &
                      (usrdata['GeneIDCount_All'] > 0)].copy()  # should keep WL's
    temp_df['LabelFLAG'] = labelix
    return temp_df

def get_gene_capacity(genes_df, database, col='GeneID'):
    """Get gene capcaity from the stored metadata"""
    capacity = (database.groupby('geneid').capacity.mean()
                .to_frame(name='GeneCapacity'))
    genes_df = genes_df.merge(capacity, how='left', left_on='GeneID', right_index=True)
    return genes_df

def get_gene_info(genes_df, database, col='GeneID'):

    subset = ['geneid', 'homologene', 'description', 'symbol', 'taxon']
    genecapacity = (database.groupby('geneid')['capacity']
                    .mean()
                    .rename('capacity_mean')
    )
    geneinfo = (database[subset]
                .drop_duplicates('geneid')
                .set_index('geneid')
                .join(genecapacity)
                .rename(columns=dict(gi = 'ProteinGI',
                                     homologene = 'HomologeneID',
                                     taxon = 'TaxonID',
                                     description = 'Description',
                                     ref = 'ProteinAccession',
                                     symbol = 'GeneSymbol',
                                     capacity_mean = 'GeneCapacity'
                ))
    )
    geneinfo.index = geneinfo.index.astype(str)
    # geneinfo['TaxonID'] = geneinfo.TaxonID.astype(str)
    out = genes_df.merge(geneinfo, how='left', left_on='GeneID', right_index=True)
    return out




def get_peptides_for_gene(genes_df, temp_df):

    full = (temp_df.groupby('GeneID')['sequence_lower']
            .agg((lambda x: frozenset(x), 'nunique'))
            .rename(columns={'<lambda>': 'PeptideSet', 'nunique': 'PeptideCount'})
            # .agg(full_op)
            .assign(PeptidePrint = lambda x: x['PeptideSet'].apply(sorted).str.join('_'))
    )
    full['PeptideSet'] = full.apply(lambda x : frozenset(x['PeptideSet']), axis=1)

    q_uniq = 'GeneIDCount_All == 1'
    q_strict = 'PSM_IDG < 4'
    q_strict_u = '{} & {}'.format(q_uniq, q_strict)

    try:
        uniq = (temp_df.query(q_uniq)
                .groupby('GeneID')['sequence_lower']
                .agg('nunique')
                .to_frame('PeptideCount_u2g'))
    except IndexError:
        uniq = pd.DataFrame(columns=['PeptideCount_u2g'])

    try:
        strict = (temp_df.query(q_strict)
                  .groupby('GeneID')['sequence_lower']
                  .agg('nunique')
                  .to_frame('PeptideCount_S'))
    except IndexError:
        strict = pd.DataFrame(columns=['PeptideCount_S'])

    try:
        s_u2g = (temp_df.query(q_strict_u)
                 .groupby('GeneID')['sequence_lower']
                 .agg('nunique')
                 .to_frame('PeptideCount_S_u2g'))
    except IndexError:
        s_u2g = pd.DataFrame(columns=['PeptideCount_S_u2g'])

    result = pd.concat((full, uniq, strict, s_u2g), copy=False, axis=1).fillna(0)
    ints = ['' + x for x in ('PeptideCount', 'PeptideCount_u2g', 'PeptideCount_S',
                                 'PeptideCount_S_u2g')]
    result[ints] = result[ints].astype(np.integer)

    genes_df = genes_df.merge(result, how='left',
                              left_on='GeneID', right_index=True)
    return genes_df

def get_psms_for_gene(genes_df, temp_df):
    psmflag = 'PSM_UseFLAG'

    total = temp_df.groupby('GeneID')[psmflag].sum()
    total.name = 'PSMs'

    q_uniq = 'GeneIDCount_All == 1'
    total_u2g = (temp_df.query(q_uniq)
                 .groupby('GeneID')[psmflag]
                 .sum())
    total_u2g.name = 'PSMs_u2g'

    q_strict = 'PSM_IDG < 4'
    total_s = (temp_df.query(q_strict)
               .groupby('GeneID')[psmflag]
               .sum())
    total_s.name = 'PSMs_S'

    q_strict_u = '{} & {}'.format(q_uniq, q_strict)
    total_s_u2g = (temp_df.query(q_strict_u)
                   .groupby('GeneID')[psmflag]
                   .sum())
    total_s_u2g.name = 'PSMs_S_u2g'

    result = (pd.concat( (total, total_u2g, total_s, total_s_u2g), copy=False, axis=1)
              .fillna(0)
              .astype(np.integer))
    genes_df = genes_df.merge(result, how='left',
                              left_on='GeneID', right_index=True)
    return genes_df


def calculate_full_areas(genes_df, temp_df, area_col, normalize):
    """ Calculates full (non distributed ) areas for gene ids.
    calculates full, gene count normalized, unique to gene,
    and unique to gene with no miscut areas.
    """
    qstring = 'AUC_UseFLAG == 1'
    full = temp_df.query(qstring).groupby('GeneID')[area_col].sum()/normalize
    full.name = 'AreaSum_max'

    # full_adj = (temp_df.query(qstring)
    #             .assign(gpAdj = lambda x: x[area_col] / x['GeneIDCount_All'])
    #             .groupby('GeneID')['gpAdj']  # temp column
    #             .sum()/normalize
    #             )
    # full_adj.name = 'AreaSum_gpcAdj'

    # qstring_s = qstring + ' & IDG < 4'
    # strict = temp_df.query(qstring_s).groupby('GeneID')[area_col].sum()
    # strict.name = ''

    qstring_u = qstring + ' & GeneIDCount_All == 1'
    uniq = temp_df.query(qstring_u).groupby('GeneID')[area_col].sum()/normalize
    uniq.name = 'AreaSum_u2g_all'

    qstring_u0 = qstring_u + ' & MissedCleavages == 0'
    uniq_0 = temp_df.query(qstring_u0).groupby('GeneID')[area_col].sum()/normalize
    uniq_0.name = 'AreaSum_u2g_0'
    result = pd.concat( (full, uniq, uniq_0), copy=False, axis=1) .fillna(0)
    genes_df = genes_df.merge(result, how='left',
                              left_on='GeneID', right_index=True)
    return genes_df

def _distribute_area(inputdata, genes_df, area_col, taxon_totals=None, taxon_redistribute=True):
    """Row based normalization of PSM area (mapped to a specific gene).
    Normalization is based on the ratio of the area of unique peptides for the
    specific gene to the sum of the areas of the unique peptides for all other genes
    that this particular peptide also maps to.
    """
    # if inputdata.AUC_UseFLAG == 0:
    #     return 0
    inputvalue = inputdata[area_col]
    geneid = inputdata['GeneID']
    gene_inputdata = genes_df.query('GeneID == @geneid')
    u2g_values = gene_inputdata['AreaSum_u2g_all'].values

    if len(u2g_values) == 1:
        u2g_area = u2g_values[0] # grab u2g info, should always be
    #of length 1
    elif len(u2g_values) > 1 :
        warn('DistArea is not singular at GeneID : {}'.format(
             datetime.now(),inputdata['GeneID']))
        distArea = 0
        # this should never happen (and never has)
    else :
        distArea = 0
        print('No distArea for GeneID : {}'.format(inputdata['GeneID']))
    # taxon_ratio = taxon_totals.get(inputdata.gene_taxon_map, 1)

    if u2g_area != 0 :
        totArea = 0
        gene_list = inputdata.GeneIDs_All.split(SEP)
        all_u2gareas = (genes_df[genes_df['GeneID'].isin(gene_list)]
                        .query('PeptideCount_u2g > 0')  # all geneids with at least 1 unique pept
                        .AreaSum_u2g_all)
        if len(all_u2gareas) > 1 and any(x == 0 for x in all_u2gareas):
            # special case with multiple u2g peptides but not all have areas, rare but does happen
            u2g_area = 0  # force to distribute by gene count (and taxon percentage if appropriate)
        else:
            totArea = all_u2gareas.sum()
            distArea = (u2g_area/totArea) * inputvalue
        #ratio of u2g peptides over total area
    elif all(gene_inputdata.IDSet == 3):
        return 0
    if u2g_area == 0:  # no uniques, normalize by genecount
        taxon_percentage = taxon_totals.get(str(inputdata.TaxonID), 1)
        distArea = inputvalue
        if taxon_percentage < 1:
            distArea *=  taxon_percentage
        gpg_selection = genes_df.GPGroup == gene_inputdata.GPGroup.values[0]
        try:
            if taxon_redistribute:
                taxonid_selection = genes_df.TaxonID == gene_inputdata.TaxonID.values[0]
                distArea /= len( genes_df[(gpg_selection) & (taxonid_selection)])
            else:
                distArea /= len( genes_df[(gpg_selection)
                ])
        except ZeroDivisionError:
            pass

    return distArea

def distribute_area(temp_df, genes_df, area_col, taxon_totals, taxon_redistribute=True):
   """Distribute psm area based on unique gene product area
   Checks for AUC_UseFLAG==1 for whether or not to use each peak for quantification
   """

   q = 'AUC_UseFLAG == 1 & GeneIDCount_All > 1'
   distarea = 'PrecursorArea_dstrAdj'
   temp_df[distarea] = 0
   # temp_df[distarea] = (temp_df.query(q)
   #                      .apply(
   #                          _distribute_area, args=(genes_df,
   #                                                      area_col,
   #                                                      taxon_totals,
   #                                                      taxon_redistribute),
   #                          axis=1)
   # )
   temp_df[distarea] = (temp_df.query(q)
                        .pipe(apply_by_multiprocessing,
                              _distribute_area,
                              workers=WORKERS,
                              func_args=(genes_df, area_col, taxon_totals, taxon_redistribute),
                              axis=1)
   )
   one_id = (temp_df.GeneIDCount_All == 1) & (temp_df.AUC_UseFLAG == 1)
   temp_df.loc[ one_id , distarea ] = temp_df.loc[ one_id, area_col ]
   temp_df[distarea].fillna(0, inplace=True)

   return

def _set2_or_3(row, genes_df, allsets):

    peptset = row.PeptideSet
    # allsets = genes_df.PeptideSet.unique()  # calculate outside this function for performance boost
    if six.PY2 and any(set(peptset) < x for x in allsets):
            return 3

    elif six.PY3 and any(peptset < allsets):
        return 3


    # check if is set 3 across multiple genes, or is set2
    gid = row.GeneID

    # sel = genes_df[ (genes_df.IDSet == 1) &
    #                 (genes_df.PeptideSet & peptset) ].query('GeneID != @gid')

    sel = genes_df[(genes_df.PeptideSet & peptset) ].query('GeneID != @gid')
    sel_idset1 = sel.query('IDSet == 1')

    in_pop = sel.PeptideSet
    in_pop_set1 = sel_idset1.PeptideSet

    in_row = sel.apply( lambda x: peptset - x['PeptideSet'], axis=1 )

    in_pop_all = set(in_pop.apply(tuple).apply(pd.Series).stack().unique())

    if not in_pop_set1.empty:
        in_pop_all_set1 = set(in_pop_set1.apply(tuple).apply(pd.Series).stack().unique())
    else:
        in_pop_all_set1 = set()

    diff = (peptset - in_pop_all)  # check if is not a subset of anything

    diff_idset1 = (peptset - in_pop_all_set1)  # check if is not a subset of set1 ids

    if len( diff_idset1 ) == 0:  # is a subset of idset1 ids
        return 3

    elif len( diff ) > 0: # is not a subset of anything
        return 2

    else:

        sel_not_idset1 = sel.query('IDSet != 1')

        if any(sel_not_idset1.PeptideSet == peptset):
            return 2  # shares all peptides with another, and is not a subset of anything

        # need to check if is a subset of everything combined, but not a subset of one thing
        # ex:
        #        PEPTIDES
        # row    =  A   B
        # match1 =  A
        # match2 =      B
        if (all( (peptset - sel.PeptideSet).apply(bool) ) and
            not all( (sel_not_idset1.PeptideSet - peptset).apply(bool) )
        ):
            return 2
        else:
            pept_lengths = sel_not_idset1.PeptideSet.apply(len)
            if len(peptset) >= pept_lengths.max():
                return 2
            else:
                return 3

            # len_shared = sel_not_idset1.PeptideSet.apply(lambda x: x & peptset).apply(len)
            # max_shared = len_shared.max()
            # all_shared_pepts = (set([x for y in sel_not_idset1.PeptideSet.values for x in y])
            #                     & peptset)


    return 3

class _DummyDataFrame:
    def eat_args(self, *args, **kwargs):
        return None
    def __getattr__(self, name):
        if name not in self.__dict__:
            return self.eat_args

def check_length_in_pipe(df):
    """Checks the length of a DataFrame in a pipe
    and if zero returns an object to suppress all errors,
    just returning None (ideally)
    """
    if len(df) == 0:
        return _DummyDataFrame()
    return df

def assign_gene_sets(genes_df, temp_df):
    all_ = genes_df.PeptideSet.unique()
    allsets = genes_df.PeptideSet.unique()
    genes_df.loc[genes_df.PeptideCount_u2g > 0, 'IDSet'] = 1
    genes_df.loc[genes_df.PeptideCount_u2g == 0, 'IDSet'] = \
                            (genes_df.query('PeptideCount_u2g == 0')
                             .pipe(check_length_in_pipe)
                             # .apply(_set2_or_3, args=(genes_df, allsets), axis=1))
                             .pipe(apply_by_multiprocessing, _set2_or_3, genes_df=genes_df, allsets=allsets,
                                   axis=1, workers=WORKERS)
                            )
    genes_df['IDSet'] = genes_df['IDSet'].fillna(3).astype(np.int8)
    # if u2g count greater than 0 then set 1
    gpg = (temp_df.groupby('GeneID')
           .PSM_IDG.min()
           .rename('IDGroup'))
    gpg_u2g = (temp_df.query('GeneIDCount_All==1')
               .groupby('GeneID')
               .PSM_IDG.min()
               .rename('IDGroup_u2g'))
    gpgs = (pd.concat([gpg, gpg_u2g], axis=1).fillna(0).astype(np.int8)
            .assign(GeneID = lambda x: x.index)
            )
    genes_df = pd.merge(genes_df, gpgs, on='GeneID', how='left')
    return genes_df

def calculate_gene_dstrarea(genes_df, temp_df, normalize=1):
    """Calculate distributed area for each gene product"""
    result = (temp_df.query('AUC_UseFLAG == 1')
              .groupby('GeneID')['PrecursorArea_dstrAdj']
              .sum()
              .divide(normalize)
              .to_frame(name='AreaSum_dstrAdj')
    )
    genes_df = genes_df.merge(result, how='left',
                              left_on='GeneID', right_index=True)
    genes_df.loc[genes_df['IDSet'] == 3, 'AreaSum_dstrAdj'] = 0
    return genes_df

def calculate_gene_razorarea(genes_df, temp_df, normalize):
    """Calculate razor area for each gene product"""

    separate_groups = lambda gpg_all : set(float(x.strip()) for z in
                                           (y.split(SEP) for y in gpg_all.values)
                                           for x in z
    )

    def razor_area(temp_df, genes_df):
        gid = temp_df.GeneID
        gpgs = (genes_df[genes_df.GeneID==gid].GPGroups_All
                .pipe(separate_groups))
        if len(gpgs) == 0:
            return 0
        allgenes = genes_df[ genes_df.GPGroup.isin(gpgs) ]
        max_uniq = allgenes.PeptideCount_u2g.max()
        gene_selection = allgenes[allgenes.PeptideCount_u2g==max_uniq]
        if len(gene_selection) != 1 or max_uniq == 0:  # no uniques
            return 0
        if gid == gene_selection.GeneID.values[0]:
            return temp_df.SequenceArea
        else:
            return 0

    temp_df['RazorArea'] = 0
    q = 'AUC_UseFLAG == 1 & GeneIDCount_All > 1'
    temp_df['RazorArea'] = (temp_df.query(q)
                                .apply(razor_area,
                                       args=(genes_df,),
                                       axis=1
                                ))
    one_id = temp_df.GeneIDCount_All == 1
    temp_df.loc[ one_id , 'RazorArea' ] = temp_df.loc[ one_id,
                                                           'SequenceArea' ]
    result = temp_df.groupby('GeneID')['RazorArea'].sum() / normalize
    result = (temp_df.groupby('GeneID')['RazorArea']
              .sum()
              .divide(normalize)
              .to_frame(name='AreaSum_razor')
    )
    genes_df = genes_df.merge(result, how='left',
                              left_on='GeneID', right_index=True)
    return genes_df

def set_gene_gpgroups(genes_df, temp_df):
    """Assign GPGroups"""

    genes_df.sort_values(by=['PSMs'], ascending=False, inplace=True)
    genes_df.reset_index(inplace=True)

    set1s =  genes_df.query('IDSet==1')

    gpg_vals = range(1, len(set1s) + 1 )
    set1_gpgs = (pd.Series(data=gpg_vals, index=set1s.index, name='GPGroup')
                 .to_frame())

    next_gpg = len(set1_gpgs) + 1
    set2_filtered = genes_df.query('IDSet==2')
    if len(set2_filtered) > 0:
    # if len(set2_filts = (set2_filtered
        set2s = (set2_filtered
                 .groupby('PeptideSet')
                 .first()
        .index)
        set2_gpgdict = (pd.Series(index=set2s,
                                  data=range(next_gpg, next_gpg + len(set2s)))
        .to_dict())
        genes_df['GPGroup'] = genes_df.PeptideSet.map(set2_gpgdict)
    else:
        genes_df['GPGroup'] = [np.nan] * len(genes_df)
    genes_df.loc[set1_gpgs.index, 'GPGroup'] = set1_gpgs


    # (temp_df[['Sequence', 'GeneIDs_All']]
    #  .drop_duplicates('Sequence')
    #  .assign(geneids = lambda x: (x['GeneIDs_All'].str.split(',')))
    # )


    gene_pept_mapping = defaultdict(set)
    pept_group_mapping = defaultdict(set)
    valid_gpgroups = genes_df[~(genes_df.GPGroup.isnull())]
    for ix, row in valid_gpgroups.iterrows():
        for pept in row['PeptideSet']:
            gene_pept_mapping[row.GeneID].add(pept)
            pept_group_mapping[pept].add(row.GPGroup)
    for ix, row in genes_df[genes_df.GPGroup.isnull()].iterrows():  # need to add set3s
        for pept in row['PeptideSet']:
            gene_pept_mapping[row.GeneID].add(pept)

    def gpg_all(gid, gene_pept_mapping, pept_group_mapping):
        gpgs = set()
        for pept in gene_pept_mapping.get(gid):
            mapping = pept_group_mapping.get(pept)
            if mapping is None:
                continue
            gpgs |= pept_group_mapping.get(pept)
        return SEP.join(str(int(x)) for x in sorted(gpgs))

    genes_df['GPGroups_All'] = genes_df.apply(lambda x: gpg_all(x['GeneID'],
                                                                    gene_pept_mapping,
                                                                    pept_group_mapping),
    axis=1)
    # do not need to multiprocess this
    # genes_df['GPGroups_All'] = genes_df.GeneID.pipe(apply_by_multiprocessing,
    #                                                 gpg_all,
    #                                                 func_args=(gene_pept_mapping, pept_group_mapping),
    #                                                 workers=WORKERS,
    # )

    genes_df['GPGroup'] = genes_df['GPGroup'].fillna(0).astype(np.int)
    # genes_df['GPGroup'].replace(to_replace='', value=float('NaN'),
    #                                 inplace=True)  # can't sort int and
    #strings, convert all strings to NaN
    return genes_df

def get_labels(usrdata, labels, labeltype='none'):
    '""labels is a dictionary of lists for each label type""'
    if labeltype == 'none':  # label free
        return ['none']
    mylabels = labels.get(labeltype)
    if mylabels is None:
        return ['none']
    included_labels = [label for label in mylabels if label in usrdata.columns]
    return included_labels

def redistribute_area_isobar(temp_df, label, labeltypes, area_col, labeltype):
    """for tmt/itraq"""
    # with_reporter = temp_df[temp_df['QuanUsage'] == 'Use']
    q = '.*{}.*'.format(labeltype.lower())
    with_reporter = temp_df[ temp_df['SequenceModi'].str.contains(q)]
    reporter_area = with_reporter[label] * with_reporter[area_col] / with_reporter[labeltypes].sum(1)
    new_area_col = '' + area_col + '_split'
    reporter_area.name = new_area_col
    temp_df = temp_df.join(reporter_area, how='right')  # only keep peptides with good reporter ion
    temp_df[new_area_col].fillna(temp_df[area_col], inplace=True)
    return temp_df, new_area_col

def concat_isobar_e2gs(rec, run, search, outdir, cols=None, labeltype=None):
    pat = re.compile('^{}_{}_{}_{}_\d+_e2g.tab'.format(rec, run, search, labeltype))
    files = list()
    for entry in os.scandir(outdir):
        if entry.is_file() and pat.search(entry.name):
            files.append(os.path.join(outdir, entry.name))
    if len(files) == 0:
        warn('No output for {}_{}_{}'.format(rec, run, search))
        return
    df = pd.concat([pd.read_table(f) for f in files])
    outf = '{}_{}_{}_{}_all_e2g.tab'.format(rec, run, search, labeltype)
    df.to_csv(os.path.join(outdir, outf), columns=cols,
                    index=False, encoding='utf-8', sep='\t')
    print('Export of {} e2g file : {}'.format(labeltype, outf))

def assign_sra(df):

    df['SRA'] = 'A'
    # cat_type = CategoricalDtype(categories=['S', 'R', 'A'],
    #                             ordered=True)
    # df['SRA'] = df['SRA'].astype(cat_type)
    df['SRA'] = df['SRA'].astype('category', categories=['S', 'R', 'A'],
                                 ordered=True)

    df.loc[ (df['IDSet'] == 1) &
            (df['IDGroup_u2g'] <= 3), 'SRA'] = 'S'

    df.loc[ (df['IDSet'] == 2) &
            (df['IDGroup'] <= 3), 'SRA'] = 'S'

    df.loc[ (df['IDSet'] == 1) &
            (df['SRA'] != 'S') &
            (df['IDGroup_u2g'] <= 5), 'SRA'] = 'R'

    df.loc[ (df['IDSet'] == 2) &
            (df['SRA'] != 'S') &
            (df['IDGroup'] <= 5), 'SRA'] = 'R'

    return df

# from ._orig_code import *
# from ._orig_code import (extract_peptideinfo, _extract_peptideinfo,
#                          _peptidome_matcher, peptidome_matcher)
def grouper(usrdata, outdir='', database=None,
            gid_ignore_file='', labels=dict(),
            contaminant_label='__CONTAMINANT__'):
    """Function to group a psm file from PD after Mascot Search"""

    def print_log_msg(df_or_None=None, msg='', *args, **kwargs):
        """Print and log a message.
        Can be used in pipe operation as always returns the
        first argument untouched.
        """
        fmt = dict(now = datetime.now(),
                   msg = msg
        )
        log_msg = '{now} | {msg}'.format(**fmt)
        print(log_msg)
        logging.info(log_msg)
        usrdata.to_logq(log_msg)
        return df_or_None

    #import RefseqInfo
    usrfile = usrdata.datafile
    # file with entry of gene ids to ignore for normalizations
    gid_ignore_list = []
    usrdata_out = usrdata.output_name('psms', ext='tab')
    if gid_ignore_file is not None and os.path.isfile(gid_ignore_file):
        print('Using gene filter file for normalization.')
        gid_ignore_list = get_gid_ignore_list(gid_ignore_file)

    area_col = 'PrecursorArea'  # set default
    normalize = 10**9

    print('Starting Grouper for exp file {}'.format(usrfile))
    print('\nFilter values set to : {}'.format(usrdata.filterstamp))
    usrdata.to_logq('{} | Starting {} for file : {}'.format(time.ctime(),
                                                            program_title, usrfile))
    # ==================== Populate gene info ================================ #
    # gene_metadata = extract_metadata(usrdata.df.metadatainfo)
    gene_taxon_dict = gene_taxon_mapper(database)
    gene_symbol_dict = gene_symbol_mapper(database)
    gene_desc_dict = gene_desc_mapper(database)
    gene_hid_dict = gene_hid_mapper(database)
    gene_protgi_dict = gene_protgi_mapper(database)
    gene_protref_dict = gene_protref_mapper(database)

    # ==================== Populate gene info ================================ #

    # usrdata.df['HID'], usrdata.df['ProteinGI'] = '', ''
    # potentially filled in later, but  likely not. Will be kept as list.

    usrdata.to_logq('{} | Finished matching PSMs to {} '\
                    'refseq'.format(time.ctime(),
                                    usrdata.taxonid))
    logging.info('Starting Grouper for exp number'\
                 '{}'.format(repr(usrdata)))
    nomatches = usrdata.df[usrdata.df['GeneIDCount_All'] == 0].Sequence  # select all
    #PSMs that didn't get a match to the refseq peptidome
    matched_psms = len(usrdata.df[usrdata.df['GeneIDCount_All'] != 0])
    unmatched_psms = len(nomatches)
    usrdata.to_logq('{} | Total identified PSMs : {}'.format(time.ctime(),
                                                             matched_psms))
    usrdata.to_logq('{} | Total unidentified PSMs : {}'.format(time.ctime(),
                                                               unmatched_psms))
    for missing_seq in nomatches:
        logging.warning('No match for sequence {} in {}'.format(missing_seq,
                                                                usrfile))
        # Store all of these sequences in the big log file, not per experiment.
    # usrdata.df = (usrdata.df.pipe(split_multiplexed)
    usrdata.df = (usrdata.df.pipe(assign_IDG, filtervalues=usrdata.filtervalues)
                  .assign(sequence_lower = lambda x: x['Sequence'].str.lower())
                  .sort_values(by=['SpectrumFile', area_col,
                                   'Sequence', 'Modifications',
                                   'Charge', 'PSM_IDG', 'IonScore', 'PEP',
                                   'q_value'],
                               ascending=[0, 0, 1, 1, 1, 1, 0, 1, 1])
                  .pipe(redundant_peaks)  # remove ambiguous peaks
                  .pipe(sum_area)
                  .pipe(auc_reflagger)  # remove duplicate sequence areas
                  .pipe(flag_AUC_PSM, usrdata.filtervalues, contaminant_label=contaminant_label,
                        phospho=usrdata.phospho)
                  .pipe(split_on_geneid)
                  .assign(TaxonID = lambda x: x['GeneID'].map(gene_taxon_dict),
                          Symbol = lambda x: x['GeneID'].map(gene_symbol_dict),
                          Description = lambda x: x['GeneID'].map(gene_desc_dict),
                          ProteinGIs = lambda x: x['GeneID'].map(gene_protgi_dict),
                          ProteinRefs = lambda x: x['GeneID'].map(gene_protref_dict),
                  )
    )
    labeltypes = get_labels(usrdata.df, labels, usrdata.labeltype)
    additional_labels = list()
    # ======================== Plugin for multiple taxons  ===================== #
    taxon_ids = usrdata.df['TaxonID'].replace(['0', 0], np.nan).dropna().unique()
    taxon_totals = dict()
    usrdata.to_logq("TaxonIDs: {}".format(len(taxon_ids)))
    # usrdata.to_logq(str(usrdata.df))
    if len(taxon_ids) == 1 or usrdata.no_taxa_redistrib:  # just 1 taxon id present
        for tid in taxon_ids: # taxon_ids is a set
            taxon_totals[tid] = 1
    elif len(taxon_ids) > 1:  # more than 1 taxon id
        taxon_totals = multi_taxon_splitter(taxon_ids, usrdata.df,
                                            gid_ignore_list, area_col)
        print('Multiple taxons found, redistributing areas...')
        usrdata.to_logq('{} | Multiple taxons found, redistributing areas.'.format(time.ctime()))
        usrdata.taxon_ratio_totals.update(taxon_totals)
        for taxon, ratio in taxon_totals.items():
            print('For the full data : {} = {}'.format(taxon, ratio))
            usrdata.to_logq('For the full data : {} = {}'.format(taxon, ratio))

    # none/SILAC loop
    # labeltypes = ['nolabel', 'heavy']  # right now only using nolabel, but in
                                       # future this can grow


    gpgcount, genecount, ibaqtot, = 0, 0, 0



    orig_area_col = 'SequenceArea'
    # Don't use the aggregated SequenceArea for TMT experiments
    if usrdata.labeltype in ('TMT', 'iTRAQ'):
        isobar_output = pd.DataFrame()  # instead of merging, we will concat
    for label in labeltypes:  # increase the range to go through more label types
        labelix = labelflag.get(label, 0)
        area_col = orig_area_col
        usrdata.to_logq('{} | Starting gene assignment for label {}.'.format(
            time.ctime(), label))
        # ==========Select only peptides flagged  with good quality=========== #
        temp_df = select_good_peptides(usrdata.df, labelix)
        if len(temp_df) == 0:  # only do if we actually have peptides selected
            print('No good peptides found for {}'.format(labelix))
            continue
        if usrdata.labeltype in ('TMT', 'iTRAQ'):
            isobar_area_col = 'PrecursorArea'  # we always use Precursor Area
            temp_df, area_col = redistribute_area_isobar(temp_df, label,
                                                         labeltypes,
                                                         area_col=isobar_area_col,
                                                         labeltype=usrdata.labeltype)
            if len(taxon_ids) > 1 and not usrdata.no_taxa_redistrib:  # more than 1 taxon id
                print('Calculating taxon ratios for label {}'.format(label))
                usrdata.to_logq('{} | Calculating taxon ratios for label{}.'.format(time.ctime(), label))
                taxon_totals = multi_taxon_splitter(taxon_ids, temp_df,
                                                    gid_ignore_list, area_col='PrecursorArea_split')
                for taxon, ratio in taxon_totals.items():
                    print('For label {} : {} = {}'.format(label, taxon, ratio))
                    usrdata.to_logq('For label {} : {} = {}'.format(label, taxon, ratio))
        elif usrdata.labeltype == 'SILAC':
            raise NotImplementedError('No support for SILAC experiments yet.')
        # ==================================================================== #

        genedata_out = usrdata.output_name(str(labelix)+'_e2g', ext='tab')
        print('{}: Populating gene table for {}.'.format(datetime.now(),
                                                         usrdata.datafile))

        msg_areas = 'Calculating peak areas for {}'.format(usrdata.datafile)
        msg_sets = 'Assigning gene sets and groups for {}'.format(usrdata.datafile)

        genes_df = (create_df(temp_df, label)
                    # .assign(TaxonID = lambda x : x['GeneID'].map(gene_taxon_dict))
                    .pipe(get_gene_info, database)
                    # .pipe(get_gene_capacity, database)
                    .pipe(get_peptides_for_gene, temp_df)
                    .pipe(get_psms_for_gene, temp_df)
                    .pipe(print_log_msg, msg=msg_areas)
                    .pipe(calculate_full_areas, temp_df, area_col, normalize)
                    .fillna(0)
                    .pipe(print_log_msg, msg=msg_sets)
                    .pipe(assign_gene_sets, temp_df)
                    .pipe(set_gene_gpgroups, temp_df)
                    .assign(Symbol = lambda x: x['GeneID'].map(gene_symbol_dict),
                            Description = lambda x: x['GeneID'].map(gene_desc_dict),
                            HIDs = lambda x: x['GeneID'].map(gene_hid_dict),
                            ProteinGIs = lambda x: x['GeneID'].map(gene_protgi_dict),
                            ProteinRefs = lambda x: x['GeneID'].map(gene_protref_dict),
                    )
                    .pipe(assign_sra)
        )


        msg_areas = 'Calculating distributed area ratio for {}'.format(usrdata.datafile)
        print_log_msg(df=None, msg=msg_areas)
        distribute_area(temp_df, genes_df, area_col, taxon_totals,
                            not usrdata.no_taxa_redistrib)
        now = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
        genes_df = (genes_df
                    .pipe(calculate_gene_dstrarea, temp_df, normalize)
                    # .pipe(calculate_gene_razorarea, temp_df, normalize)
                    .assign(iBAQ_dstrAdj = lambda x:
                            np.divide(x['AreaSum_dstrAdj'],
                                      x['GeneCapacity']))
                    .sort_values(by=['GPGroup'], ascending=True)
                    .reset_index()
                    .assign(EXPRecNo = usrdata.recno,
                            EXPRunNo = usrdata.runno,
                            EXPSearchNo = usrdata.searchno,
                            AddedBy = usrdata.added_by,
                            CreationTS = now,
                            ModificationTS = now)
                    .sort_values(by=['IDSet', 'GPGroup'])
        )

        # =============================================================#

        genecount += len(genes_df)
        ibaqtot += genes_df[~genes_df.GeneID.isin(
            gid_ignore_list)].iBAQ_dstrAdj.sum()
        if usrdata.labeltype in ('TMT', 'iTRAQ'):
            isobar_output = pd.concat([isobar_output, temp_df])
        genes_df.to_csv(os.path.join(usrdata.outdir, genedata_out), columns=E2G_COLS,
                        index=False, encoding='utf-8', sep='\t')
        usrdata.to_logq('{} | Export of genetable for labeltype {}'\
                        'completed.'.format(
                            time.ctime(),
                            label))
    # ========================================================================= #
    usrdata.df.drop('metadatainfo', axis=1, inplace=True)  # Don't need this
                                                           # column anymore.
    if len(usrdata.df) == 0:
        print('No protein information for {}.\n'.format(repr(usrdata)))
        usrdata.to_logq('No protein information for {}.\n'.format(repr(usrdata)))
        usrdata.flush_log()
        return

    usrdata.df['ModificationTS'] = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
    #usrdata.df['HIDs'] = ''  # will be populated later
    #usrdata.df['HIDCount'] = ''  # will be populated later


    usrdata.to_logq('{} | Starting peptide ranking.'.format(time.ctime()))
    dstr_area = 'PrecursorArea_dstrAdj'
    area_col_to_use = dstr_area if dstr_area in usrdata.df.columns else orig_area_col
    data_cols = DATA_COLS
    if usrdata.labeltype in ('TMT', 'iTRAQ'):
        if usrdata.labeltype == 'TMT':
            data_cols = DATA_COLS + ['TMT_126', 'TMT_127_N', 'TMT_127_C', 'TMT_128_N',
                                     'TMT_128_C', 'TMT_129_N', 'TMT_129_C', 'TMT_130_N',
                                     'TMT_130_C', 'TMT_131', 'QuanInfo', 'QuanUsage']
        elif usrdata.labeltype == 'iTRAQ':
            data_cols = DATA_COLS + ['iTRAQ_114', 'iTRAQ_115', 'iTRAQ_116', 'iTRAQ_117',
                                     'QuanInfo', 'QuanUsage']

        isobar_output.reset_index(inplace=True)
        rank_df = pd.DataFrame()
        if 'LabelFLAG' not in isobar_output:
            isobar_output['LabelFLAG'] = np.nan
            isobar_output['PeptRank'] = 0
        for label in isobar_output.LabelFLAG.dropna().unique():
            q = 'LabelFLAG == @label'
            isobar_rank = rank_peptides(isobar_output.query(q),
                                        area_col=area_col_to_use,
                                        ranks_only=True)
            rank_df = pd.concat([rank_df, isobar_rank.to_frame()])
        now = datetime.now().strftime("%m/%d/%Y %H:%M:%S")

        isobar_output = (isobar_output.join(rank_df, how='left')
                         .assign(PeptRank = lambda x: (x['PeptRank']
                                                           .fillna(0)
                                                           .astype(np.integer)),
                                 ModificationTS = now,
                                 HIDs = '',  # will be populated later
                                 HIDCount = ''  # will be populated later
                      )
        )
    else:
        usrdata.df = (pd.merge(usrdata.df, temp_df, how='left')
                      .pipe(rank_peptides, area_col=area_col_to_use)
                      .assign(PrecursorArea_split = lambda x: x['PrecursorArea'],
                              PeptRank = lambda x: (x['PeptRank']
                                                        .fillna(0)
                                                        .astype(np.integer))
                      )
        )
        # usrdata.df['PrecursorArea_split'] = usrdata.df['PrecursorArea']
                              # didn't get a rank gets a rank of 0
    msg = 'Peptide ranking complete for {}.'.format(usrdata.datafile)
    print_log_msg(msg=msg)

    export_metadata(program_title=program_title, usrdata=usrdata,
                    matched_psms=matched_psms, unmatched_psms=unmatched_psms,
                    usrfile=usrfile, taxon_totals=usrdata.taxon_ratio_totals,
                    outname=usrdata.output_name('metadata', ext='json'), outpath=usrdata.outdir)

    msfname = usrdata.output_name('msf', ext='tab')
    msfdata = spectra_summary(usrdata)
    msfdata.to_csv(os.path.join(usrdata.outdir, msfname), index=False, sep='\t')

    if usrdata.labeltype in ('TMT', 'iTRAQ'):
        if not all(x in isobar_output.columns.values for x in set(data_cols) - set(_EXTRA_COLS)):
            print('Potential error, not all columns filled.')
            print([x for x in data_cols if x not in isobar_output.columns.values])
        data_cols = [x for x in data_cols if x in isobar_output.columns.values]
        concat_isobar_e2gs(usrdata.recno, usrdata.runno, usrdata.searchno,
                           usrdata.outdir, cols=E2G_COLS, labeltype=usrdata.labeltype)
        # usrdata.df = pd.merge(usrdata.df, temp_df, how='left')
        isobar_output.to_csv(os.path.join(usrdata.outdir, usrdata_out), columns=data_cols,
                             index=False, encoding='utf-8', sep='\t')
    else:
        if not all(x in usrdata.df.columns.values for x in set(data_cols) - set(_EXTRA_COLS)):
            print('Potential error, not all columns filled.')
            print([x for x in data_cols if x not in usrdata.df.columns.values])
        data_cols = [x for x in data_cols if x in usrdata.df.columns.values]
        usrdata.df.to_csv(os.path.join(usrdata.outdir, usrdata_out), columns=data_cols,
                          index=False, encoding='utf-8', sep='\t')

    usrdata.to_logq('{} | Export of datatable completed.'.format(time.ctime())+
                    '\nSuccessful grouping of file completed.')
    usrdata.flush_log()

    print('Successful grouping of {} completed.\n'.format(repr(usrdata)))

def calculate_breakup_size(row_number):
    return ceil(row_number/4)

def set_modifications(usrdata):

    to_replace = {'DeStreak' : 'des', 'Deamidated' : 'dam', 'Carbamidomethyl' : 'car',
                  'Oxidation' : 'oxi', 'Phospho' : 'pho',
                  'Acetyl': 'ace', 'GlyGly' : 'gg', 'Label:13C(6)' : 'lab',
                  'Label:13C(6)+GlyGly' : 'labgg',
                  '\)\(': ':'}
    modis_abbrev = usrdata.Modifications.replace(regex=to_replace).fillna('')
    modis_abbrev.name = 'Modifications_abbrev'
    usrdata = usrdata.join(modis_abbrev)
    modifications = usrdata.apply(lambda x :
                                  seq_modi(x['Sequence'],
                                           x['Modifications_abbrev'],
                                           to_replace.values()
                                  ),
                                  axis=1
    )
    (usrdata['Sequence'], usrdata['SequenceModi'],
     usrdata['SequenceModiCount'], usrdata['LabelFLAG']) = zip(*modifications)
    return usrdata

def load_fasta(refseq_file):
    REQUIRED_COLS = ('geneid', 'sequence')
    ADDITIONAL_COLS = ('description', 'gi', 'homologene', 'ref', 'taxon', 'symbol')
    gen = fasta_dict_from_file(refseq_file)

    # routine for converting things to integer if possible, else leave as string
    l = []
    for g in gen:
        if 'geneid' not in g:
            continue
        if g['geneid'] == '' or (isinstance(g['geneid'], (int, float)) and g['geneid'].isnull()):
            continue
        for k,v in g.items():
            if g[k] == 'nan':
                g[k] = -1
            try:
                g[k] = int(v)
                g['']
            except Exception as e:
                pass
        l.append(g)
    # end
    df = (pd.DataFrame.from_dict(l)
          # .replace(np.nan, '')
    )
    if not all(x in df.columns for x in REQUIRED_COLS):
        missing = ', '.join(x for x in REQUIRED_COLS if x not in df.columns)
        fmt = 'Invalid FASTA file : {} is missing the following identifiers : {}\n'
        raise ValueError(fmt.format(refseq_file, missing))
    missing_cols = (x for x in ADDITIONAL_COLS if x not in df.columns)
    for col in missing_cols:
        df[col] = ''
    return df

ENZYME = {'trypsin': dict(cutsites=('K', 'R'), exceptions=None),
          'trypsin/P': dict(cutsites=('K', 'R'), exceptions=('P',)),
}
# TODO: add the rest
 # 'trypsin/P', 'chymotrypsin', 'LysC', 'LysN', 'GluC', 'ArgC' 'AspN',}

def _match(usrdatas, refseq_file, miscuts=2, enzyme='trypsin/P'):

    enzyme_rule = ENZYME[enzyme]
    print(u'Using peptidome {} with rule {}'.format(refseq_file, enzyme))

    # database = pd.read_table(refseq_file, dtype=str)
    # rename_refseq_cols(database, refseq_file)
    database = load_fasta(refseq_file)
    database['capacity'] = 1
    breakup_size = calculate_breakup_size(len(database))
    counter = 0
    prot = defaultdict(list)
    for ix, row in database.iterrows():
        counter += 1
        fragments, fraglen = protease(row.sequence, minlen=7,
                                      # cutsites=['K', 'R'],
                                      # exceptions=['P'],
                                      miscuts=miscuts,
                                      **enzyme_rule
        )
        database.loc[ix, 'capacity'] = fraglen
        for fragment in fragments: # store location in the DataFrame for the peptide's parent
            prot[fragment].append(ix)

        if counter > breakup_size:
            for usrdata in usrdatas:
                usrdata.df = peptidome_matcher(usrdata.df, prot)  # match peptides to peptidome
            counter = 0
            del prot # frees up memory, can get quite large otherwise
            prot = defaultdict(list)
    else:
        for usrdata in usrdatas:
            usrdata.df  = peptidome_matcher(usrdata.df, prot)  # match peptides to peptidome
        del prot

    # now extract info based on index
    for usrdata in usrdatas:
        if usrdata.searchdb is None:
            usrdata.searchdb = refseq_file
        extract_peptideinfo(usrdata, database)

    return database

def match(usrdatas, refseqs, enzyme='trypsin/P'):
    """
    Match psms with fasta database
    Input is list of UserData objects and an optional dictionary of refseqs
    """
    inputdata_refseqs = set([usrdata.taxonid for usrdata in usrdatas])
    databases = dict()
    sortfunc = lambda x: x.taxon_miscut_id
    usrdata_sorted = sorted(usrdatas, key=sortfunc)
    for k, g in itertools.groupby(usrdata_sorted, key=sortfunc):
        group = list(g)
        taxonid = group[0].taxonid
        miscuts = group[0].miscuts
        refseq = refseqs.get(taxonid)

        if refseq is None:
            err = 'No refseq file available for {}'.format(taxonid)
            warn(err)
            for u in group:
                u.ERROR = err
                u.EXIT_CODE = 1
            continue

        database = _match(group, refseq, miscuts=miscuts, enzyme=enzyme)
        databases[taxonid] = database
    # for organism in refseqs:
    #     if any(x == int(organism) for x in inputdata_refseqs):
    #                                     # check if we even need to read the
    #                                     # peptidome for that organism
    #         database = _match([usrdata for usrdata in usrdatas if usrdata.taxonid==organism],
    #                           refseqs[organism])
    #         databases[organism] = database

    return databases

def column_identifier(df, aliases):
    column_names = dict()
    for col in aliases:
        for alias in aliases[col]:
            name = [dfcolumn for dfcolumn in df.columns if dfcolumn==alias]
            if len(name)==1:
                column_names[col] = name[0]
                break
    return column_names

# Idea
# REQUIRED_HEADERS = ['Sequence', 'Modifications', 'PrecursorArea',
#                     'Charge', 'IonScore', 'q_value', 'PEP', 'SpectrumFile',
#                     'RTmin', 'DeltaMassPPM']
# def check_required_headers(df):
#     if not all(x in df.columns for x in REQUIRED_HEADERS):
#         missing = [x for x in REQUIRED_HEADERS if x not in df.columns]
#         fmt = 'Invalid input file, missing {}'.format(', '.join(missing))
#         raise ValueError(fmt)


def set_up(usrdatas, column_aliases):
    """Set up the usrdata class for analysis
    Read data, rename columns (if appropriate), populate base data"""
    for usrdata in usrdatas:
        EXIT_CODE = usrdata.read_csv(sep='\t')  # read from the stored psms file
        if EXIT_CODE != 0:
            print('Error with reading {!r}'.format(usrdata))
            continue
        if column_aliases:
            standard_names = column_identifier(usrdata.df, column_aliases)
            protected_names = ('Modified sequence')
            usrdata.df.rename(columns={v: k
                                       for k,v in standard_names.items()},
                              inplace=True
            )
            redundant_cols = [x for x in usrdata.df.columns if
                              (x not in standard_names.keys() and
                               x not in protected_names)
            ]
            # print(usrdata.df.memory_usage().sum())
            usrdata.df = usrdata.df.drop(redundant_cols, axis=1)
            # print(usrdata.df.memory_usage().sum())

        # usrdata.df = usrdata.populate_base_data()
        usrdata.populate_base_data()
        if 'DeltaMassPPM' not in usrdata.df:
            usrdata.df['DeltaMassPPM'] = 0
        if 'SpectrumFile' not in usrdata.df:
            usrdata.df['SpectrumFile'] = None
        if 'RTmin' not in usrdata.df:
            usrdata.df['RTmin'] = 0
        if 'PEP' not in usrdata.df.columns:
            usrdata.df['PEP'] = 0  # not trivial to make a category due to sorting
            # usrdata.categorical_assign('PEP', 0, ordered=True)
        # check_required_headers(usrdata.df)


        if 'MissedCleavages' not in usrdata.df.columns:
            usrdata.df['MissedCleavages'] =\
                                        usrdata.df.apply(lambda x:\
                                                         calculate_miscuts(x['Sequence'],
                                                                           targets=('K', 'R')),
                                                         axis=1)
        if not 'q_value' in usrdata.df.columns:
            try:
                usrdata.df['q_value'] = usrdata.df['PEP'] / 10  # rough approximation
            except KeyError:
                warn('No q_value available')
                usrdata.df['q_value'] = 0
        if not 'PSMAmbiguity' in usrdata.df.columns:
            # usrdata.df['PSMAmbiguity'] = 'Unambiguous'
            usrdata.categorical_assign('PSMAmbiguity', 'Umambiguous')  # okay?
        if not usrdata.pipeline == 'MQ':  # MaxQuant already has modifications
            usrdata.df = set_modifications(usrdata.df)
        else:
            # usrdata.df['SequenceModi'] = usrdata.df['Modified sequence']
            usrdata.df.rename(columns={'Modified sequence': 'SequenceModi'},
                              inplace=True)
            usrdata.df['SequenceModiCount'] = count_modis_maxquant(usrdata.df)
            usrdata.categorical_assign('LabelFlag', 0) #TODO: handle this properly
            # usrdata.df['LabelFLAG'] = 0  #TODO: handle this properly
    return

def rename_refseq_cols(df, filename):
    fasta_h = ['TaxonID', 'HomologeneID', 'GeneID', 'ProteinGI', 'FASTA']
    # regxs = {x: re.compile(x) for x in _fasta_h}
    to_rename = dict()
    for header in fasta_h:
        matched_col = [col for col in df.columns
                       if header in col]
        if len(matched_col) != 1:
            raise ValueError("""Could not identify the correct
            column in the reference sequence database for '{}'
            for file : {}""".format(header, filename))
        to_rename[matched_col[0]] = header
    df.rename(columns=to_rename, inplace=True)

WORKERS = 1

def main(usrdatas=[], fullpeptread=False, inputdir='', outputdir='', refs=dict(),
         rawfilepath=None, column_aliases=dict(), gid_ignore_file='', labels=dict(),
         raise_on_error=False, contaminant_label='__CONTAMINANT__', enzyme='trypsin/P', workers=1):
    """
    refs :: dict of taxonIDs -> refseq file names
    """
    # ====================Configuration Setup / Loading======================= #
    global WORKERS
    WORKERS = workers
    if imagetitle:
        fancyprint(program_title, 12)  # ascii rt
        #  fancyprint('Malovannaya lab',10)  #
    elif not imagetitle:
        print(program_title)  # not as exciting as ascii art

    print('\nrelease date: {}'.format(__copyright__))
    print('Python version ' + sys.version)
    print('Pandas version: ' + pd.__version__)

    startTime = datetime.now()

    print('\nStart at {}'.format(startTime))
    # logging.info('Start at {}'.format(startTime))

    # first set the modifications. Importantly fills in X with the predicted amino acid
    set_up(usrdatas, column_aliases)
    # for ud in usrdatas:
    #     print(ud, ud.EXIT_CODE)
    if all(usrdata.EXIT_CODE != 0 for usrdata in usrdatas):
        return usrdatas

    databases = match([x for x in usrdatas if x.EXIT_CODE == 0], refs, enzyme=enzyme)
    if all(usrdata.EXIT_CODE != 0 for usrdata in usrdatas):
        return usrdatas

    # failed_exps = []
    for usrdata in usrdatas:
        if usrdata.EXIT_CODE != 0:
            continue
        try:
            grouper(usrdata,
                    database=databases[usrdata.taxonid],
                    gid_ignore_file=gid_ignore_file, labels=labels,
                    contaminant_label=contaminant_label)
            usrdata.EXIT_CODE = 0
        except Exception as e:  # catch and store all exceptions, won't crash
                                # the whole program at least
            usrdata.EXIT_CODE = 1
            usrdata.ERROR = traceback.format_exc()
            stack = traceback.format_exc()
            # failed_exps.append((usrdata, e))
            usrdata.to_logq('Failure for file of experiment {}.\n'\
                  'The reason is : {}'.format(repr(usrdata), e))
            usrdata.to_logq(stack)
            print(stack)
            # for s in stack:
            #     usrdata.to_logq(s, sep='')
            #     print(s, sep='')
            usrdata.flush_log()
            print('Failure for file of experiment {}.\n'\
                  'The reason is : {}'.format(repr(usrdata), e))

            logging.warn('Failure for file of experiment {}.\n'\
                         'The reason is : {}'.format(repr(usrdata), e))
            if raise_on_error:
                raise  # usually don't need to raise, will kill the script. Re-enable
                       # if need to debug and find where errors are
    print('Time taken : {}\n'.format(datetime.now() - startTime))
    return usrdatas
