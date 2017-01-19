#===============================================================================#
# PyGrouper - Alex Saltzman
import re, os, sys, time
import json
import logging
from collections import defaultdict
from functools import partial
from math import ceil
from warnings import warn
from configparser import ConfigParser
from itertools import repeat

import numpy as np
import pandas as pd

from . import _version
from .subfuncts import *

from ._orig_code import timed
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
}

try:
    from PIL import Image, ImageFont, ImageDraw
    imagetitle = True
except ImportError:
    imagetitle = False

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
    PSM_count    = len(data[data.psm_PSM_UseFLAG==1])
    dmass_median = data.DeltaMassPPM.median()
    return(RT_min, RT_max, IonScore_min, IonScore_max, q_min, q_max,
           PEP_min, PEP_max, area_min, area_max, PSM_count, dmass_median)

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
    msfdata['RawFileName']       = list(set(usrdata.df.SpectrumFile.tolist()))
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
     msfdata['DeltaMassPPM_med']) = list(zip(*summary_info))

    rawfile_info = msfdata.apply(lambda x:
                                    _get_rawfile_info(usrdata.rawfiledir,
                                                      x['RawFileName']),
                                 axis=1)

    msfdata['RawFileSize'], msfdata['RawFileTS'] = list(zip(*rawfile_info))

    msfdata.rename(columns={c: 'msf_'+c for c in msfdata.columns}, inplace=True)

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

def _extract_peptideinfo(row):

    # row = database.loc[ix]
    if len(row) == 0:
        return ('', 0, '', 0, '', 0, '', 0, ())
    result = (
        # ','.join(row['GeneID'].dropna().unique()),
        ','.join(str(x) for x in set(row['GeneID'])),
        row['GeneID'].nunique(dropna=True),

        # ','.join(row['TaxonID'].dropna().unique()),
        ','.join(str(x) for x in set(row['TaxonID'])),
        row['TaxonID'].nunique(dropna=True),

        # ','.join(row['ProteinGI'].dropna().unique()),
        ','.join(str(x) for x in set(row['ProteinGI'])),
        row['ProteinGI'].nunique(dropna=True),

        # ','.join(row['HomologeneID'].dropna().unique()),
        ','.join(str(x) for x in set(row['HomologeneID'])),
        row['HomologeneID'].nunique(dropna=True),
        tuple(row['capacity']),
          # row['capacity'].mean(),
    )
    return result



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
    # t0= time.time()

    info = ixs.apply(lambda x : _extract_peptideinfo(database.loc[x]))
    (usrdata.df['psm_GeneList'],
     usrdata.df['psm_GeneCount'],
     usrdata.df['psm_TaxonIDList'],
     usrdata.df['psm_TaxonCount'],
     usrdata.df['psm_ProteinList'],
     usrdata.df['psm_ProteinCount'],
     usrdata.df['psm_HIDList'],
     usrdata.df['psm_HIDCount'],
     usrdata.df['psm_GeneCapacities']) = zip(*info)
    # print(time.time()-t0, '\n')
    usrdata.df['psm_TaxonIDList'] = usrdata.df['psm_TaxonIDList'].dropna().astype(str)
    usrdata.df['psm_HIDList'] = usrdata.df['psm_HIDList'].fillna('')

    return 0

def gene_taxon_mapper(df):
    """Returns a dictionary with mapping:
    gene -> taxon
    Input is the metadata extracted previously"""
    return {x[1].GeneID: x[1].TaxonID for x in df.iterrows()}



def assign_IDG(usrdata):
    ion_score_bins = usrdata.filtervalues.get('ion_score_bins', (10, 20, 30))
    df = usrdata.df
    df['psm_PSM_IDG'] = pd.cut(df['IonScore'],
                               bins=(0, *ion_score_bins, np.inf),
                               labels=[7, 5, 3, 1], include_lowest=True,
                               right=False).astype('int')
    df.loc[ df['q_value'] > .01, 'psm_PSM_IDG' ] += 1
    return usrdata

def make_seqlower(usrdata, col='Sequence'):
    """Make a new column called sequence_lower from a DataFrame"""
    usrdata['sequence_lower'] = usrdata[col].str.lower()
    return usrdata

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
            drop_duplicates(subset=['SpectrumFile','psm_SequenceModi', 'Charge', 'PrecursorArea'])
    peaks.is_copy = False  # duplicate dataframe in memory
    peaks['psm_Peak_UseFLAG'] = 1
    usrdata = usrdata.join(peaks['psm_Peak_UseFLAG'])
    usrdata['psm_Peak_UseFLAG'] = usrdata.psm_Peak_UseFLAG.fillna(0)
    print('Redundant peak areas removed : ', len(usrdata)-len(peaks))
    return usrdata

def sum_area(usrdata, area_col):
    """Sum the area of similar peaks
    New column psm_SequenceArea is created
    """
    usrdata['Sequence_set'] = usrdata['Sequence'].apply(lambda x: tuple(set(list(x))))
    summed_area = pd.DataFrame(usrdata[usrdata.psm_Peak_UseFLAG==1][
        ['psm_SequenceModi', 'Charge',
         'PrecursorArea']].groupby(['psm_SequenceModi', 'Charge'])["PrecursorArea"].sum())
    # summed_area = usrdata[usrdata.psm_Peak_UseFLAG==1][area_col]

    summed_area.reset_index(inplace=True)
    summed_area.rename(columns={area_col: 'psm_SequenceArea'}, inplace=True)
    usrdata = usrdata.merge(summed_area, how='left', on=['psm_SequenceModi', 'Charge'])
    #usrdata = usrdata.join(summed_area['psm_SequenceArea'], how='left')
    return usrdata

def auc_reflagger(usrdata, area_col):
    """Remove duplicate sequence areas
    area_col input is original area column,
    area_col output is renamed to psm_SequenceArea
    """
    #usrdata['Sequence_set'] = usrdata['Sequence'].apply(lambda x: tuple(set(list(x))))
    no_dups = usrdata.sort_values(by=['psm_SequenceModi', 'Charge', 'psm_SequenceArea',
                                      'psm_PSM_IDG', 'IonScore', 'PEP', 'q_value'],
                                  ascending=[1,1,0,1,0,1,1]).drop_duplicates(subset=
                                                                   ['psm_SequenceArea',
                                                                    'Charge',
                                                                    'psm_SequenceModi',])
    no_dups.is_copy = False
    no_dups['AUC_reflagger'] = 1
    usrdata = usrdata.join(no_dups[['AUC_reflagger']])
    usrdata['AUC_reflagger'] = usrdata['AUC_reflagger'].fillna(0)
    # area_col = 'psm_SequenceArea'
    return usrdata

def export_metadata(program_title='version',usrdata=None, matched_psms=0, unmatched_psms=0,
                    usrfile='file', taxon_totals=dict(), outname=None, outpath='.', **kwargs):
    """Update iSPEC database with some metadata information
    """
    print('{} | Updating experiment runs table in iSPEC.'.format(time.ctime()))
    #print('Number of matched psms : ', matched_psms)
    d = dict(
        version=program_title,
        searchdb=usrdata.searchdb,
        filterstamp=usrdata.filterstamp,
        matched_psms=matched_psms,
        unmatched_psms=unmatched_psms,
        inputname=usrdata.datafile,
        hu=taxon_totals.get(9606, 0),
        mou=taxon_totals.get(10090, 0),
        gg=taxon_totals.get(9031, 0),
        recno=usrdata.recno,
        runno=usrdata.runno,
        searchno=usrdata.searchno
    )
    with open(os.path.join(outpath, outname), 'w') as f:
        json.dump(d, f)

def split_on_geneid(usrdata):
    """Duplicate psms based on geneids. Areas of each psm is recalculated based on
unique peptides unique for its particular geneid later.
"""
    glstsplitter = usrdata['psm_GeneList'].str.split(',').apply(pd.Series,
                                                                1).stack()
    glstsplitter.name = 'psm_GeneID'  # give the series a name
    glstsplitter = pd.DataFrame(glstsplitter)
    # use multi index to determine the original row
    glstsplitter['psm_oriFLAG'] = glstsplitter.index.map(lambda x: 1 if x[-1] == 0 else 0)

    glstsplitter.index = glstsplitter.index.droplevel(-1)  # get rid of
                                                           # multi-index

    usrdata = usrdata.join(glstsplitter)  # usrdata gains column 'psm_GeneID'
                                          #from glstsplitter Series
    usrdata.reset_index(inplace=True, drop=True)  # drop=True ?
    return usrdata

def rank_peptides(usrdata, area_col, ranks_only=False):
    """Rank peptides here
    area_col is sequence area_calculator
    ranks_only returns just the ranks column. This does not reset the original index
    """

    usrdata = usrdata.sort_values(by=['psm_GeneID', area_col,
                                      'psm_SequenceModi',
                                      'Charge', 'psm_PSM_IDG', 'IonScore', 'PEP',
                                      'q_value'],
                                  ascending=[1, 0, 0, 1, 1, 0, 1, 1])
    if not ranks_only:  # don't reset index for just the ranks
        usrdata.reset_index(inplace=True)  # drop=True ?
    usrdata.Modifications.fillna('', inplace=True)  # must do this to compare nans
    usrdata[area_col].fillna(0, inplace=True)  # must do this to compare
    #nans

    grouped = usrdata[ (usrdata.psm_AUC_UseFLAG == 1) & \
                       (usrdata.psm_PSM_UseFLAG == 1) & \
                       (usrdata.psm_Peak_UseFLAG == 1)    ].groupby(['psm_GeneID', 'psm_LabelFLAG'])
    ranks = grouped.cumcount() + 1  # add 1 to start the peptide rank at 1, not 0
    ranks.name = 'psm_PeptRank'
    if ranks_only:
        return ranks

    usrdata = usrdata.join(ranks)

    return usrdata


def flag_AUC_PSM(df, fv):
    if fv['pep'] =='all' : fv['pep'] = float('inf')
    if fv['idg'] =='all' : fv['idg'] = float('inf')
    df['psm_AUC_UseFLAG'] = 1
    df['psm_PSM_UseFLAG'] = 1
    df.loc[(df['Charge'] < fv['zmin']) | (df['Charge'] > fv['zmax']),
           ['psm_AUC_UseFLAG', 'psm_PSM_UseFLAG']] = 0

    df.loc[df['psm_SequenceModiCount'] > fv['modi'],
           ['psm_AUC_UseFLAG', 'psm_PSM_UseFLAG']] = 0

    df.loc[(df['IonScore'].isnull() & df['q_value'].isnull()),
           ['psm_AUC_UseFLAG', 'psm_PSM_UseFLAG']] = 1, 0

    df.loc[df['IonScore'] < fv['ion_score'],
           ['psm_AUC_UseFLAG', 'psm_PSM_UseFLAG']] = 0

    df.loc[df['q_value'] > fv['qvalue'],
           ['psm_AUC_UseFLAG', 'psm_PSM_UseFLAG']] = 0

    df.loc[df['PEP'] > fv['pep'],
           ['psm_AUC_UseFLAG', 'psm_PSM_UseFLAG']] = 0
    df.loc[df['psm_PSM_IDG'] > fv['idg'],
           ['psm_AUC_UseFLAG', 'psm_PSM_UseFLAG']] = 0

    df.loc[(df['psm_Peak_UseFLAG'] == 0) & (df['PSMAmbiguity'].str.lower()=='unambiguous'),
           ['psm_AUC_UseFLAG', 'psm_PSM_UseFLAG']] = 1
    df.loc[(df['psm_Peak_UseFLAG'] == 0) & (df['PSMAmbiguity'].str.lower()!='unambiguous'),
           ['psm_AUC_UseFLAG', 'psm_PSM_UseFLAG']] = 0

    df.loc[ df['AUC_reflagger'] == 0, 'psm_AUC_UseFLAG'] = 0
    return df


def gene_taxon_map(usrdata, gene_taxon_dict):
    """make 'gene_taxon_map' column per row which displays taxon for given gene"""
    usrdata['psm_TaxonID'] = usrdata['psm_GeneID'].map(gene_taxon_dict)
    return usrdata


def get_all_taxons(taxonidlist):
    """Return a set of all taxonids from
    usrdata.psm_TaxonIDList"""

    taxon_ids = set(','.join(x for x in taxonidlist
                             if x.strip()).split(','))
    return taxon_ids

def multi_taxon_splitter(taxon_ids, usrdata, gid_ignore_list, area_col):
    """Plugin for multiple taxons
    Returns a dictionary with the totals for each detected taxon"""
    taxon_totals = dict()
    for taxon in taxon_ids:
        taxon = str(taxon)
        #all_others = [x for x in taxon_ids if x != taxon]
        uniq_taxon = usrdata[
            #(usrdata._data_tTaxonIDList.str.contains(taxon)) &
            #(~usrdata._data_tTaxonIDList.str.contains('|'.join(all_others)))&
            (usrdata['psm_TaxonIDList'] == str(taxon)) &
            #(usrdata['psm_PSM_IDG']<9) &  # this is redunant with AUC_UseFLAG
            (~usrdata['psm_GeneID'].isin(gid_ignore_list)) &
            (usrdata['psm_AUC_UseFLAG'] == 1)
        ]
        taxon_totals[taxon] = (uniq_taxon[area_col] / uniq_taxon['psm_GeneCount']).sum()
        tot_unique = sum(taxon_totals.values())  #sum of unique
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

def create_e2g_df(inputdf, label, inputcol='psm_GeneID'):
    """Create and return a DataFrame with gene/protein information from the input
    peptide DataFrame"""
    return pd.DataFrame({'e2g_GeneID':
                         list(set(inputdf[inputcol])),
                         'e2g_EXPLabelFLAG': labelflag.get(label, 0)})

def select_good_peptides(usrdata, labelix):
    """Selects peptides of a given label with the correct flag and at least one genecount"""
    temp_df = usrdata[# (usrdata['psm_AUC_UseFLAG'] == 1) &
                      (usrdata['psm_PSM_UseFLAG'] == 1) &
                      (usrdata['psm_GeneCount'] > 0)].copy()  # should keep WL's
    temp_df['psm_LabelFLAG'] = labelix
    return temp_df

def get_gene_capacity(genes_df, database, col='e2g_GeneID'):
    """Get gene capcaity from the stored metadata"""
    capacity = database.groupby('GeneID').capacity.mean()
    capacity.name = 'e2g_GeneCapacity'
    genes_df = genes_df.merge(capacity.to_frame(), how='left', left_on='e2g_GeneID', right_index=True)
    return genes_df


def get_peptides_for_gene(genes_df, temp_df):
    full_op = {'e2g_PeptideSet': 'unique', 'e2g_PeptideCount': 'nunique'}
    strict_op = {'e2g_PeptideCount_S': 'nunique'}
    u2g_op = {'e2g_PeptideCount_u2g': 'nunique'}
    s_u2g_op = {'e2g_PeptideCount_S_u2g': 'nunique'}
    full = (temp_df.groupby('psm_GeneID')['sequence_lower']
            .agg(full_op)
            .assign(e2g_PeptidePrint = lambda x: x['e2g_PeptideSet'].str.join('_'))
    )
    full['e2g_PeptideSet'] = full.apply(lambda x : set(x['e2g_PeptideSet']), axis=1)

    q_uniq = 'psm_GeneCount == 1'
    q_strict = 'psm_PSM_IDG < 4'
    q_strict_u = '{} & {}'.format(q_uniq, q_strict)

    uniq = (temp_df.query(q_uniq)
            .groupby('psm_GeneID')['sequence_lower']
            .agg(u2g_op))

    strict = (temp_df.query(q_strict)
              .groupby('psm_GeneID')['sequence_lower']
              .agg(strict_op))
    s_u2g= (temp_df.query(q_strict_u)
            .groupby('psm_GeneID')['sequence_lower']
            .agg(s_u2g_op))

    result = pd.concat( (full, uniq, strict, s_u2g), copy=False, axis=1) .fillna(0)
    genes_df = genes_df.merge(result, how='left',
                              left_on='e2g_GeneID', right_index=True)
    return genes_df



def get_psms_for_gene(genes_df, temp_df):
    psmflag = 'psm_PSM_UseFLAG'

    total = temp_df.groupby('psm_GeneID')[psmflag].sum()
    total.name = 'e2g_PSMs'

    q_uniq = 'psm_GeneCount == 1'
    total_u2g = (temp_df.query(q_uniq)
                 .groupby('psm_GeneID')[psmflag]
                 .sum()
    )
    total_u2g.name = 'e2g_PSMs_u2g'

    q_strict = 'psm_PSM_IDG < 4'
    total_s = (temp_df.query(q_strict)
               .groupby('psm_GeneID')[psmflag]
               .sum()
    )
    total_s.name = 'e2g_PSMs_S'

    q_strict_u = '{} & {}'.format(q_uniq, q_strict)
    total_s_u2g = (temp_df.query(q_strict_u)
                   .groupby('psm_GeneID')[psmflag]
                   .sum()
    )
    total_s_u2g.name = 'e2g_PSMs_S_u2g'

    result = pd.concat( (total, total_u2g, total_s, total_s_u2g), copy=False, axis=1) .fillna(0)
    genes_df = genes_df.merge(result, how='left',
                              left_on='e2g_GeneID', right_index=True)
    return genes_df


def calculate_protein_area(genes_df, temp_df, area_col, normalize):

    qstring = 'psm_AUC_UseFLAG == 1'
    full = temp_df.query(qstring).groupby('psm_GeneID')[area_col].sum()/normalize
    full.name = 'e2g_nGPArea_Sum_max'

    full_adj = (temp_df.query(qstring)
                .assign(gpAdj = lambda x: x[area_col] / x['psm_GeneCount'])
                .groupby('psm_GeneID')['gpAdj']  # temp column
                .sum()/normalize
                )
    full_adj.name = 'e2g_nGPArea_Sum_cgpAdj'
    # qstring_s = qstring + ' & psm_PSM_IDG < 4'
    # strict = temp_df.query(qstring_s).groupby('psm_GeneID')[area_col].sum()
    # strict.name = ''

    qstring_u = qstring + ' & psm_GeneCount == 1'
    uniq = temp_df.query(qstring_u).groupby('psm_GeneID')[area_col].sum()/normalize
    uniq.name = 'e2g_nGPArea_Sum_u2g'

    qstring_u0 = qstring_u + ' & MissedCleavages == 0'
    uniq_0 = temp_df.query(qstring_u0).groupby('psm_GeneID')[area_col].sum()/normalize
    uniq_0.name = 'e2g_nGPArea_Sum_u2g_all'
    result = pd.concat( (full, full_adj, uniq, uniq_0), copy=False, axis=1) .fillna(0)
    genes_df = genes_df.merge(result, how='left',
                              left_on='e2g_GeneID', right_index=True)
    return genes_df

def _distribute_psm_area(inputdata, genes_df, area_col, taxon_totals=None, taxon_redistribute=True):
    """Row based normalization of PSM area (mapped to a specific gene).
    Normalization is based on the ratio of the area of unique peptides for the
    specific gene to the sum of the areas of the unique peptides for all other genes
    that this particular peptide also maps to.
    """

    if inputdata.psm_AUC_UseFLAG == 0:
        return 0
    inputvalue = inputdata[area_col]
    gene_inputdata = genes_df[ genes_df['e2g_GeneID'] == inputdata['psm_GeneID']]
    u2gPept = (genes_df[genes_df['e2g_GeneID']==inputdata['psm_GeneID']]
               ['e2g_nGPArea_Sum_u2g_all']).values

    if len(u2gPept) == 1: u2gPept = u2gPept[0] # grab u2g info, should always be
    #of length 1
    elif len(u2gPept) > 1 :
        warn('DistArea is not singular at GeneID : {}'.format(
             datetime.now(),inputdata['psm_GeneID']))
        distArea = 0
        # this should never happen (and never has)
    else :
        distArea = 0
        print('No distArea for GeneID : {}'.format(inputdata['psm_GeneID']))
    # taxon_ratio = taxon_totals.get(inputdata.gene_taxon_map, 1)
    if u2gPept != 0 :
        totArea = 0
        gene_list = inputdata.psm_GeneList.split(',')
        totArea = genes_df[
             genes_df['e2g_GeneID'].isin(gene_list)
                              ].e2g_nGPArea_Sum_u2g_all.sum()
        distArea = (u2gPept/totArea) * inputvalue
        #ratio of u2g peptides over total area
    elif all(gene_inputdata.e2g_IDSet == 3):
        distArea = 0
    elif u2gPept == 0:  # no uniques, normalize by genecount
        taxon_percentage = taxon_totals.get(inputdata.psm_TaxonID, 1)
        distArea = inputvalue
        if taxon_percentage < 1:
            distArea *=  taxon_percentage
        try: # still needs work
            if taxon_redistribute:
                distArea /= len( genes_df[(genes_df.e2g_GPGroup == gene_inputdata.e2g_GPGroup.values[0]) &
                                          (genes_df.e2g_TaxonID == gene_inputdata.e2g_TaxonID.values[0])
                ])
            else:
                distArea /= len( genes_df[(genes_df.e2g_GPGroup == gene_inputdata.e2g_GPGroup.values[0])
                ])
            # distArea = inputvalue/inputdata.psm_GeneCount
        except ZeroDivisionError:
            pass
        # distArea = inputvalue/inputdata.psm_GeneCount

    return distArea

def distribute_psm_area(temp_df, genes_df, area_col, taxon_totals, taxon_redistribute=True):
   """Distribute psm area based on unique gene product area"""

   temp_df['psm_PrecursorArea_dstrAdj'] = temp_df.apply(
       _distribute_psm_area, args=(genes_df,
                                   area_col,
                                   taxon_totals,
                                   taxon_redistribute),
   axis=1)

   return temp_df


def _assign_gene_sets(genes_df,genes_df_all, temp_df ):
    IDquery, Pquery = genes_df[['e2g_GeneID','e2g_PeptideSet']]

    if any(temp_df[temp_df.psm_GeneID==IDquery]['psm_GeneCount']==1):
        idset = 1
    elif not any(genes_df_all.e2g_PeptideSet.values > Pquery):
        idset = 2
    elif any(genes_df_all.e2g_PeptideSet.values > Pquery):
        idset = 3
    try:
        idgroup = min(temp_df[temp_df.psm_GeneID == IDquery].psm_PSM_IDG)
    except ValueError:
        idgroup = 0
    try :
        idgroup_u2g  = min(temp_df[(temp_df.psm_GeneID == IDquery) & (temp_df.psm_GeneCount == 1)].psm_PSM_IDG)
    except ValueError :
        idgroup_u2g = 0

    return idset, idgroup, idgroup_u2g

def assign_gene_sets(genes_df, temp_df):
    """Assign IDSet and IDGroup"""
    genesets = genes_df.apply(_assign_gene_sets, args=(genes_df, temp_df,), axis=1)
    (genes_df['e2g_IDSet'], genes_df['e2g_IDGroup'],
     genes_df['e2g_IDGroup_u2g']) = list(zip(*genesets))

    return genes_df

# def __assign_gene_sets(genes_df, temp_df):
#     import ipdb; ipdb.set_trace()
#     set1s = temp_df[temp_df['psm_GeneCount']==1].psm_GeneID
#     genes_df['e2g_IDSet'] = 3
#     genes_df['e2g_IDGroup_u2g'] = 0

#     genes_df.loc[ genes_df.e2g_GeneID.isin(set1s), 'e2g_IDSet'] = 1
#     set_counts = (genes_df['e2g_PeptideSet']
#                   .value_counts()
#     )
#     no_subsets = set_counts[set_counts == 1]


#     return genes_df




def calculate_gene_dstrarea(genes_df, temp_df, normalize):
    """Calculate distributed area for each gene product"""
    g = temp_df.groupby('psm_GeneID')
    result = g['psm_PrecursorArea_dstrAdj'].sum() / normalize
    result.name = 'e2g_nGPArea_Sum_dstrAdj'
    genes_df = genes_df.merge(result.to_frame(), how='left',
                              left_on='e2g_GeneID', right_index=True)
    genes_df.loc[genes_df['e2g_IDSet'] == 3, 'e2g_nGPArea_Sum_dstrAdj'] = 0
    return genes_df


def _GPG_helper(idset,peptideset, df_all,last):
    if idset == 3:
        gpg =  ''
    else:
        if idset == 1:
            gpg =  last + 1
        elif idset == 2:
            selection = df_all[df_all.e2g_PeptideSet.values == peptideset].e2g_GPGroup
            selection = selection[selection != '']
            if len(selection) == 1:
                gpg =  selection.values[0]
            elif len(selection) == 0:
                gpg =  last + 1
            elif len(selection) > 1:
                uniq_sel = []
                for k in range(len(selection)):
                    if selection.iloc[k] not in uniq_sel:
                        uniq_sel.append(selection.iloc[k])
                    if len(uniq_sel) == 1:
                        gpg =  uniq_sel[0]
                    elif len(uniq_sel) > 1:
                        warn('More than one IDSet type 2 genes '
                              'already has assignment')
                        gpg =  genes_df['e2g_GPGroup']

    return gpg, True if gpg < last else False

def _GPG_all_helper(genes_df,df_all):

    GPGall = set()
    for pept in genes_df.e2g_PeptideSet:
        shared_values = [value for value in
                         df_all[df_all['e2g_PeptidePrint'].str.contains(
                              pept, regex=False, case=False)].e2g_GPGroup.values]
        # regex=False since we don't need regex, and is faster.
        for s in shared_values:
            if s !='':
                GPGall.add(s)

    return str([k for k in GPGall]).strip('[').strip(']')


def set_gene_gpgroups(genes_df):
    """Assign GPGroups"""
    genes_df['e2g_GPGroup'] = ''
    genes_df.sort_values(by=['e2g_PSMs'], ascending=False, inplace=True)
    genes_df.index = list(range(0, len(genes_df)))
    last = 0
    for i in range(len(genes_df)):  # The logic behind it makes
        #sense,
        #but the implementation is weird.
        # print(last)  # for debugging
        if genes_df.loc[i]['e2g_IDSet'] != 3:
            genes_df.loc[i,'e2g_GPGroup'], lessthan = _GPG_helper(genes_df.at[i,'e2g_IDSet'],
                                                                  genes_df.at[i,'e2g_PeptideSet'],
                                                                  genes_df, last
            )

        if isinstance(genes_df.loc[i]['e2g_GPGroup'],int) and not lessthan:
            last = genes_df.loc[i, 'e2g_GPGroup']

    genes_df['e2g_GPGroups_All'] = genes_df.apply(_GPG_all_helper,
                                                  args=(genes_df,),
                                                  axis=1)
    genes_df['e2g_GPGroup'].replace(to_replace='', value=float('NaN'),
                                    inplace=True)  # can't sort int and
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

def redistribute_area_tmt(temp_df, label, labeltypes, area_col):
    """for tmt"""
    # with_reporter = temp_df[temp_df['QuanUsage'] == 'Use']
    with_reporter = temp_df[ temp_df['psm_SequenceModi'].str.contains('.*tmt.*')]
    reporter_area = with_reporter[label] * with_reporter[area_col] / with_reporter[labeltypes].sum(1)
    new_area_col = 'psm_' + area_col + '_split'
    reporter_area.name = new_area_col
    temp_df = temp_df.join(reporter_area, how='right')  # only keep peptides with good reporter ion
    temp_df[new_area_col].fillna(temp_df[area_col], inplace=True)
    return temp_df, new_area_col

def concat_tmt_e2gs(rec, run, search, outdir, cols=None):
    pat = re.compile('^{}_{}_{}_TMT_\d+_e2g.tab'.format(rec, run, search))
    files = list()
    for entry in os.scandir(outdir):
        if entry.is_file() and pat.search(entry.name):
            files.append(os.path.join(outdir, entry.name))
    if len(files) == 0:
        warn('No output for {}_{}_{}'.format(rec, run, search))
        return
    df = pd.concat([pd.read_table(f) for f in files])
    outf = '{}_{}_{}_TMT_all_e2g.tab'.format(rec, run, search)
    df.to_csv(os.path.join(outdir, outf), columns=cols,
                    index=False, encoding='utf-8', sep='\t')
    print('Export of TMT e2g file : {}'.format(outf))

# from ._orig_code import *
# from ._orig_code import (extract_peptideinfo, _extract_peptideinfo,
#                          _peptidome_matcher, peptidome_matcher)
def grouper(usrdata, outdir='', database=None,
            gid_ignore_file='', labels=dict()):
    """Function to group a psm file from PD after Mascot Search"""
    #import RefseqInfo
    usrfile = usrdata.datafile
    # file with entry of gene ids to ignore for normalizations
    #gid_ignore_file = 'pygrouper_geneignore.txt'
    gid_ignore_list = []

    usrdata_out = usrdata.output_name('psms', ext='tab')
    if gid_ignore_file is not None and os.path.isfile(gid_ignore_file):
        print('Using gene filter file for normalization.')
        gid_ignore_list = get_gid_ignore_list(gid_ignore_file)

    if usrdata.quant_source.strip() == 'AUC':
        area_col = 'PrecursorArea'  # we always use Precursor Area
        normalize = 10**9
    elif usrdata.quant_source.strip() == 'Intensity':
        area_col = 'Intensity'
        normalize = 10**5

    print('Starting Grouper for exp file {}'.format(usrfile))
    print('\nFilter values set to : {}'.format(usrdata.filterstamp))
    logfilestr = usrdata.output_name(ext='log')
    logfile = open(os.path.join(usrdata.outdir, logfilestr),
                   'w+')  # append to previously opened log file
    logfile.write('{} | Starting {} for file : {}\n'.format(
        time.ctime(),program_title, usrfile))
    # ==================== Populate gene info ================================ #
    # gene_metadata = extract_metadata(usrdata.df.metadatainfo)
    gene_taxon_dict = gene_taxon_mapper(database)

    # usrdata.df = extract_genelist(usrdata.df)


    # ==================== Populate gene info ================================ #

    usrdata.df['psm_HID'], usrdata.df['psm_ProteinGI'] = '', ''
    # potentially filled in later,
    #fields will exist in database at least

    logfile.write('{} | Finished matching PSMs to {} '\
                  'refseq\n'.format(time.ctime(),
                                    usrdata.taxonid))
    logging.info('Starting Grouper for exp number'\
                 '{}'.format(repr(usrdata)))
    nomatches = usrdata.df[usrdata.df['psm_GeneCount'] == 0].Sequence  # select all
    #PSMs that didn't get a match to the refseq peptidome
    matched_psms = len(usrdata.df[usrdata.df['psm_GeneCount'] != 0])
    unmatched_psms = len(nomatches)
    logfile.write('{} | Total identified PSMs : {}\n'.format(time.ctime(),
                                                             matched_psms))
    logfile.write('{} | Total unidentified PSMs : {}\n'.format(time.ctime(),
                                                               unmatched_psms))
    for missing_seq in nomatches:
        logging.warning('No match for sequence {} in {}'.format(missing_seq,
                                                                usrfile))
        # Store all of these sequences in the big log file, not per experiment.
    logfile.write('{} | Starting grouper.\n'.format(time.ctime()))
    usrdata = assign_IDG(usrdata)
    usrdata.df = make_seqlower(usrdata.df)
    usrdata.df = usrdata.df.sort_values(by=['SpectrumFile', area_col,
                                      'Sequence', 'Modifications',
                                      'Charge', 'psm_PSM_IDG', 'IonScore', 'PEP',
                                      'q_value'], ascending=[0, 0, 1, 1, 1, 1, 0, 1, 1])
    usrdata.df = redundant_peaks(usrdata.df)
    # remove ambiguous peaks
    usrdata.df = sum_area(usrdata.df, area_col) # area_col is PrecursorArea / Intensity
    # now remove duplicate sequence areas
    usrdata.df = auc_reflagger(usrdata.df, area_col)
    # area_col = 'psm_SequenceArea'  # this is equal to what is returned by auc_reflagger
    #usrdata.to_csv(os.path.join(outdir, usrdata_out),index=False, sep='\t')
    #print('exiting...')
    #return
    # ============= Gather all genes that all peptides map to =============== #
    usrdata.df = split_on_geneid(usrdata.df)
    # ========================================================================= #

    # Now calculate AUC and PSM use flags
    usrdata.df = flag_AUC_PSM(usrdata.df, usrdata.filtervalues)

    usrdata.df = gene_taxon_map(usrdata.df, gene_taxon_dict)
    # Flag good quality peptides
        # ======================== Plugin for multiple taxons  ===================== #
    taxon_ids = usrdata.df['psm_TaxonID'].dropna().unique()
    # taxon_ids = get_all_taxons(usrdata.df['psm_TaxonIDList'].tolist())
    # taxon_ids = get_all_taxons(usrdata.df['psm_TaxonIDList'].str.strip().dropna().tolist())
    #area_col_new = 'psm_Area_taxonAdj'
    taxon_totals = dict()
    # print(taxon_ids)
    if len(taxon_ids) == 1 or usrdata.no_taxa_redistrib:  # just 1 taxon id present
        taxon_totals[list(taxon_ids)[0]] = 1  # taxon_ids is a set
    #    usrdata[area_col_new] = usrdata[area_col]
    elif len(taxon_ids) > 1:  # more than 1 taxon id
        taxon_totals = multi_taxon_splitter(taxon_ids, usrdata.df, gid_ignore_list, area_col)
        print('Multiple taxons found, redistributing areas...')
        logfile.write('{} | Multiple taxons found, '\
                      'redistributing areas.\n'.format(time.ctime()))
    pd.options.mode.chained_assignment = None  # default='warn'

    # none/SILAC loop
    # labeltypes = ['nolabel', 'heavy']  # right now only using nolabel, but in
                                       # future this can grow
    gpgcount, genecount, ibaqtot, = 0, 0, 0
    e2g_cols = ['e2g_EXPRecNo', 'e2g_EXPRunNo', 'e2g_EXPSearchNo',
                'e2g_EXPLabelFLAG', 'e2g_AddedBy',
                'e2g_CreationTS', 'e2g_ModificationTS', 'e2g_TaxonID',
                'e2g_GeneID',
                'e2g_IDSet', 'e2g_IDGroup', 'e2g_IDGroup_u2g',
                'e2g_GPGroup', 'e2g_GPGroups_All', 'e2g_PSMs',
                'e2g_PSMs_u2g', 'e2g_PeptidePrint', 'e2g_PeptideCount',
                'e2g_PeptideCount_u2g', 'e2g_PeptideCount_S',
                'e2g_PeptideCount_S_u2g', 'e2g_nGPArea_Sum_cgpAdj',
                'e2g_nGPArea_Sum_u2g', 'e2g_nGPArea_Sum_u2g_all',
                'e2g_nGPArea_Sum_max', 'e2g_nGPArea_Sum_dstrAdj',
                'e2g_GeneCapacity', 'e2g_n_iBAQ_dstrAdj']  # cols of interest

    labeltypes = get_labels(usrdata.df, labels, usrdata.labeltype)
    additional_labels = list()

    # orig_area_col = 'psm_SequenceArea' if usrdata.labeltype == 'None' else area_col
    orig_area_col = 'psm_SequenceArea'
    # Don't use the aggregated SequenceArea for TMT experiments
    if usrdata.labeltype == 'TMT':
        psm_tmtoutput = pd.DataFrame()  # instead of merging, we will concat
    for label in labeltypes:  # increase the range to go through more label types
        labelix = labelflag.get(label, 0)
        area_col = orig_area_col
        logfile.write('{} | Starting gene assignment for label {}.\n'.format(
            time.ctime(), label))
        # ==========Select only peptides flagged  with good quality=========== #
        mylabelix = labelix
        # if usrdata.labeltype == 'TMT':
        #     mylabelix = 0 # originally would be the label type as indicated
        #                   # by the modification, but Proteome Discoverer changed how they
        #                   # mark the label.
        #                   # Really no need for marking the label (but maybe for SILAC)
        temp_df = select_good_peptides(usrdata.df, mylabelix)
        if usrdata.labeltype == 'TMT':
            if usrdata.quant_source.strip() == 'AUC':
                tmt_area_col = 'PrecursorArea'  # we always use Precursor Area
            elif usrdata.quant_source.strip() == 'Intensity':
                tmt_area_col = 'Intensity'
            temp_df, area_col = redistribute_area_tmt(temp_df, label, labeltypes, area_col=tmt_area_col)
        elif usrdata.labeltype == 'SILAC':
            raise NotImplementedError('No support for SILAC experiments yet.')
        # ==================================================================== #
        if len(temp_df) == 0:  # only do if we actually have peptides selected
            continue

        genedata_out = usrdata.output_name(labelix, 'e2g', ext='tab')
        print('{}: Populating gene table for {}.'.format(datetime.now(),
                                                            usrdata.datafile))
        # genes_df['_e2g_GeneID'] = Set(temp_df['_data_GeneID']) #only in 2.7
        genes_df = create_e2g_df(temp_df, label)

        genes_df['e2g_TaxonID'] = genes_df.apply(lambda x: gene_taxon_dict.get(x['e2g_GeneID']),
                                                 axis=1)

        genes_df = get_gene_capacity(genes_df, database)
        genes_df = get_peptides_for_gene(genes_df, temp_df)
        genes_df = get_psms_for_gene(genes_df, temp_df)

        print('{}: Calculating peak areas for {}.'.format(
            datetime.now(), usrdata.datafile))
        logging.info('{}: Calculating peak areas for {}.'.format(
            datetime.now(), usrdata.datafile))
        logfile.write('{} | Calculating peak areas.\n'.format(time.ctime()))

        genes_df = calculate_protein_area(genes_df, temp_df, area_col, normalize).fillna(0)
        # pandas may give a warning from this though it is fine
        print('{}: Calculating distributed area ratio for {}.'.format(
            datetime.now(), usrdata.datafile))
        logging.info('{}: Calculating distributed area ratio for {}.'.format(
            datetime.now(), usrdata.datafile))
        logfile.write('{} | Calculating distributed area ratio.\n'.format(
            time.ctime()))

        print('{}: Assigning gene sets and groups for {}.'.format(
            datetime.now(), usrfile))
        logging.info('{}: Assigning gene sets and groups for {}.'.format(
            datetime.now(), usrfile))
        logfile.write('{} | Assigning gene sets and groups.\n'.format(
            time.ctime()))

        genes_df = assign_gene_sets(genes_df, temp_df)

        genes_df = set_gene_gpgroups(genes_df)


        temp_df = distribute_psm_area(temp_df, genes_df, area_col, taxon_totals,
                                      not usrdata.no_taxa_redistrib)


        #genes_df['e2g_GeneCapacity'] = genes_df.e2g_GeneCapacity.astype('float')
        #genes_df._e2g_GeneCapacity.dtype)  # for debugging
        genes_df = calculate_gene_dstrarea(genes_df, temp_df, normalize)
        genes_df['e2g_n_iBAQ_dstrAdj'] = \
            genes_df.e2g_nGPArea_Sum_dstrAdj / genes_df.e2g_GeneCapacity

        genes_df.sort_values(by=['e2g_GPGroup'], ascending=True, inplace=True)
        genes_df.index = list(range(0, len(genes_df)))  # reset the index
        gpgcount += genes_df.e2g_GPGroup.max()  # do this before filling na
        genes_df['e2g_GPGroup'].fillna('', inplace=True)  # convert all NaN
                                                        #back to empty string
        # =============================================================#
        genes_df['e2g_EXPRecNo'] = usrdata.recno
        genes_df['e2g_EXPRunNo'] = usrdata.runno
        genes_df['e2g_EXPSearchNo'] = usrdata.searchno
        genes_df['e2g_EXPTechRepNo'] = usrdata.techrepno
        genes_df['e2g_AddedBy'] = usrdata.added_by
        genes_df['e2g_CreationTS'] = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
        genes_df['e2g_ModificationTS'] = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
        #renamed_genecols = [exp_setup.get(genecol, genecol) for genecol in e2g_cols]
        #torename = {k:v for k, v in exp_setup.items() if k.startswith('e2g_')}
        # get a log of gene info
        genecount += len(genes_df)
        ibaqtot += genes_df[~genes_df.e2g_GeneID.isin(
            gid_ignore_list)].e2g_n_iBAQ_dstrAdj.sum()

        #genes_df.rename(columns=torename, inplace=True)
        # print(os.path.join(outdir, genedata_out))
        if usrdata.labeltype == 'TMT':
            psm_tmtoutput = pd.concat([psm_tmtoutput, temp_df])
        genes_df.to_csv(os.path.join(usrdata.outdir, genedata_out), columns=e2g_cols,
                        index=False, encoding='utf-8', sep='\t')
        logfile.write('{} | Export of genetable for labeltype {}'\
                        'completed.\n'.format(
                            time.ctime(),
                            label))

            # ========================================================================= #

    # ----------------End of none/silac loop--------------------------------- #

    usrdata.df.drop('metadatainfo', axis=1, inplace=True)  # Don't need this
                                      # column anymore.
    if len(usrdata.df) == 0:
        print('No protein information for {}.\n'.format(repr(usrdata)))
        logfile.write('No protein information for {}.\n'.format(repr(usrdata)))
        logfile.close()
        return


    usrdata.df['psm_ModificationTS'] = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
    #usrdata.df['psm_HIDList'] = ''  # will be populated later
    #usrdata.df['psm_HIDCount'] = ''  # will be populated later
    data_cols = ['psm_EXPRecNo', 'psm_EXPRunNo', 'psm_EXPSearchNo',
                 'psm_EXPTechRepNo', 'Sequence',
                 'PSMAmbiguity', 'Modifications', 'ActivationType',
                 'DeltaScore', 'DeltaCn', 'Rank', 'SearchEngineRank',
                 'PrecursorArea', 'q_value', 'PEP',
                 'Decoy Peptides Matched',
                 'IonScore',
                 'Peptides Matched', 'MissedCleavages',
                 'IsolationInterference', 'IonInjectTime',
                 'Intensity', 'Charge', 'mzDa', 'MHDa',
                 'DeltaMassDa', 'DeltaMassPPM', 'RTmin',
                 'FirstScan', 'LastScan', 'MSOrder', 'MatchedIons',
                 'SpectrumFile', 'psm_AddedBy',
                 'psm_oriFLAG',
                 'psm_CreationTS', 'psm_ModificationTS', 'psm_GeneID',
                 'psm_GeneList', 'psm_GeneCount', 'psm_ProteinGI',
                 'psm_ProteinList', 'psm_ProteinCount',
                 'psm_HID', 'psm_HIDList', 'psm_HIDCount',
                 'psm_TaxonID', 'psm_TaxonIDList', 'psm_TaxonCount',
                 'psm_PSM_IDG', 'psm_SequenceModi',
                 'psm_SequenceModiCount', 'psm_LabelFLAG',
                 'psm_PeptRank', 'psm_AUC_UseFLAG', 'psm_PSM_UseFLAG',
                 'psm_Peak_UseFLAG', 'psm_SequenceArea', 'psm_PrecursorArea_split',
                 'psm_PrecursorArea_dstrAdj']


    logfile.write('{} | Starting peptide ranking.\n'.format(time.ctime()))
    dstr_area = 'psm_PrecursorArea_dstrAdj'
    area_col_to_use = dstr_area if dstr_area in usrdata.df.columns else orig_area_col
    if usrdata.labeltype == 'TMT':
        data_cols = data_cols + ['TMT_126', 'TMT_127_N', 'TMT_127_C', 'TMT_128_N',
                                 'TMT_128_C', 'TMT_129_N', 'TMT_129_C', 'TMT_130_N',
                                 'TMT_130_C', 'TMT_131', 'QuanInfo', 'QuanUsage']
        psm_tmtoutput.reset_index(inplace=True)
        rank_df = pd.DataFrame()
        for label in psm_tmtoutput.psm_LabelFLAG.unique():
            tmtrank = rank_peptides(psm_tmtoutput[psm_tmtoutput.psm_LabelFLAG == label],
                                    area_col=area_col_to_use, ranks_only=True)
            rank_df = pd.concat([rank_df, tmtrank])
        rank_df.columns = ['psm_PeptRank']
            # rank_df = rank_df.join(tmtrank, how='outer')
        psm_tmtoutput = psm_tmtoutput.join(rank_df, how='left')
        psm_tmtoutput['psm_PeptRank'] = psm_tmtoutput['psm_PeptRank'].fillna(0)  # anyone who
                              # didn't get a rank gets a rank of 0
    #print('Length of usrdata after merge : ',len(usrdata))

        psm_tmtoutput['psm_ModificationTS'] = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
        psm_tmtoutput['psm_HIDList'] = ''  # will be populated later
        psm_tmtoutput['psm_HIDCount'] = ''  # will be populated later
    else:
        usrdata.df = pd.merge(usrdata.df, temp_df, how='left')
        # rare case where no PSMs pass into the temp_df of quantified PSMs
        usrdata.df = rank_peptides(usrdata.df, area_col=area_col_to_use)
        usrdata.df['psm_PeptRank'] = usrdata.df['psm_PeptRank'].fillna(0)  # anyone who
        # TODO : Make this valid for AUC and Intensity based quantification
        usrdata.df['psm_PrecursorArea_split'] = usrdata.df['PrecursorArea']
                              # didn't get a rank gets a rank of 0
    #print('Length of usrdata after merge : ',len(usrdata))

    logfile.write('{} | Peptide ranking complete.\n'.format(time.ctime()))
    print('{}: Peptide ranking complete for {}.'.format(datetime.now(), usrdata.datafile))
    logging.info('{}: Peptide ranking complete for {}.'.format(datetime.now(),
                                                               usrdata.datafile))
    #usrdata.to_csv(usrdata_out, columns=usrdata.columns,
                                #encoding='utf-8', sep='\t')


    export_metadata(program_title=program_title, usrdata=usrdata, matched_psms=matched_psms,
                    unmatched_psms=unmatched_psms, usrfile=usrfile, taxon_totals=taxon_totals,
                    outname=usrdata.output_name('metadata', ext='json'), outpath=usrdata.outdir)

    msfname = usrdata.output_name('msf', ext='tab')
    msfdata = spectra_summary(usrdata)
    msfdata.to_csv(os.path.join(usrdata.outdir, msfname), index=False, sep='\t')

    if usrdata.labeltype == 'TMT':
        if not all(x in psm_tmtoutput.columns.values for x in data_cols):
            print('Potential error, not all columns filled.')
            print([x for x in data_cols if x not in psm_tmtoutput.columns.values])
            data_cols = [x for x in data_cols if x in psm_tmtoutput.columns.values]
        concat_tmt_e2gs(usrdata.recno, usrdata.runno, usrdata.searchno,
                        usrdata.outdir, cols=e2g_cols)
        usrdata.df = pd.merge(usrdata.df, temp_df, how='left')
        psm_tmtoutput.to_csv(os.path.join(usrdata.outdir, usrdata_out), columns=data_cols,
                          index=False, encoding='utf-8', sep='\t')
    else:
        if not all(x in usrdata.df.columns.values for x in data_cols):
            print('Potential error, not all columns filled.')
            print([x for x in data_cols if x not in usrdata.df.columns.values])
            data_cols = [x for x in data_cols if x in usrdata.df.columns.values]
        usrdata.df.to_csv(os.path.join(usrdata.outdir, usrdata_out), columns=data_cols,
                          index=False, encoding='utf-8', sep='\t')


    logfile.write('{} | Export of datatable completed.\n'.format(time.ctime()))
    logfile.write('Successful grouping of file completed.')
    logfile.close()

    print('Successful grouping of {} completed.\n' \
          .format(repr(usrdata)))

def calculate_breakup_size(row_number):
    return ceil(row_number/4)

def set_modifications(usrdata):

    to_replace = {'DeStreak' : 'des', 'Deamidated' : 'dam', 'Carbamidomethyl' : 'car',
                  'Oxidation' : 'oxi', 'Phospho' : 'pho',
                  'Acetyl': 'ace', 'GlyGly' : 'gg', 'Label:13C(6)' : 'lab', 'Label:13C(6)+GlyGly' : 'labgg',
                  '\)\(': ':'}
    modis_abbrev = usrdata.Modifications.replace(regex=to_replace)
    modis_abbrev.name = 'Modifications_abbrev'
    usrdata = usrdata.join(modis_abbrev)
    modifications = usrdata.apply(lambda x :
                                  seq_modi(x['Sequence'],
                                           x['Modifications_abbrev'],
                                           to_replace.values()
                                  ),
                                  axis=1
    )
    (usrdata['Sequence'], usrdata['psm_SequenceModi'],
     usrdata['psm_SequenceModiCount'], usrdata['psm_LabelFLAG']) = list(zip(*modifications))
    return usrdata

def _match(usrdatas, refseq_file):
    print('Using peptidome {} '.format(refseq_file))
    database = pd.read_table(refseq_file, dtype=str)
    rename_refseq_cols(database, refseq_file)
    database['capacity'] = 1
    breakup_size = calculate_breakup_size(len(database))
    counter = 0
    prot = defaultdict(list)
    for ix, row in database.iterrows():
        counter += 1
        fragments, fraglen = protease(row.FASTA, minlen=7,
                                      cutsites=['K', 'R'],
                                      exceptions=['P'])
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

def match(usrdatas, refseqs):
    """
    Match psms with fasta database
    Input is list of UserData objects and an optional dictionary of refseqs
    """
    inputdata_refseqs = set([usrdata.taxonid for usrdata in usrdatas])
    databases = dict()
    for organism in refseqs:
        if any(x == int(organism) for x in inputdata_refseqs):
                                        # check if we even need to read the
                                        # peptidome for that organism
            database = _match([usrdata for usrdata in usrdatas if usrdata.taxonid==organism],
                              refseqs[organism])
            databases[organism] = database

    return usrdatas, databases

def column_identifier(df, aliases):
    column_names = dict()
    for col in aliases:
        for alias in aliases[col]:
            name = [dfcolumn for dfcolumn in df.columns if dfcolumn==alias]
            if len(name)==1:
                column_names[col] = name[0]
                break
    return column_names

def set_up(usrdatas, column_aliases):
    """Set up the usrdata class for analysis

    Read data, rename columns (if appropriate), populate base data"""
    for usrdata in usrdatas:
        usrdata.read_csv(sep='\t')  # read from the stored psms file
        if column_aliases:
            standard_names = column_identifier(usrdata.df, column_aliases)
            usrdata.df.rename(columns={v: k
                                       for k,v in standard_names.items()},
                              inplace=True
            )
        if 'PrecursorArea' not in usrdata.df.columns and usrdata.quant_source == 'AUC' and 'Intensity' in usrdata.df.columns:
            # explicitly rename as MaxQuant referrs to the PSM area as Intensity
            usrdata.df.rename(columns={'Intensity': 'PrecursorArea'}, inplace=True)
        # usrdata.df = usrdata.populate_base_data()
        usrdata.populate_base_data()
        if 'MissedCleavages' not in usrdata.df.columns:
            usrdata.df['MissedCleavages'] =\
                                        usrdata.df.apply(lambda x:\
                                                         calculate_miscuts(x['Sequence'],
                                                                           targets=('K', 'R')),
                                                         axis=1)
        if not 'q_value' in usrdata.df.columns:
            usrdata.df['q_value'] = usrdata.df['PEP'] / 10  # rough approximation
        if not 'PSMAmbiguity' in usrdata.df.columns:
            usrdata.df['PSMAmbiguity'] = 'Unambiguous'
        if not usrdata.pipeline == 'MQ':  # MaxQuant already has modifications
            usrdata.df = set_modifications(usrdata.df)
        else:
            # usrdata.df['psm_SequenceModi'] = usrdata.df['Modified sequence']
            usrdata.df.rename(columns={'Modified sequence': 'psm_SequenceModi'},
                              inplace=True)
            usrdata.df['psm_SequenceModiCount'] = count_modis_maxquant(usrdata.df)
            usrdata.df['psm_LabelFLAG'] = 0  #TODO: handle this properly
    return usrdatas

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


def main(usrdatas=[], fullpeptread=False, inputdir='', outputdir='', refs=dict(),
         rawfilepath=None, column_aliases=dict(), gid_ignore_file='', labels=dict()):
    """
    refs :: dict of taxonIDs -> refseq file names
    """
    # ====================Configuration Setup / Loading======================= #
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
    usrdatas = set_up(usrdatas, column_aliases)

    usrdatas, databases = match(usrdatas, refs)

    failed_exps = []
    for usrdata in usrdatas:
        try:
            grouper(usrdata,
                    database=databases[usrdata.taxonid],
                    gid_ignore_file=gid_ignore_file, labels=labels)
        except Exception as e:  # catch and store all exceptions, won't crash
                                # the whole program at least
            failed_exps.append((usrdata, e))
            print('Failure for file of experiment {}.\n'\
                  'The reason is : {}'.format(repr(usrdata), e))
            logging.warn('Failure for file of experiment {}.\n'\
                         'The reason is : {}'.format(repr(usrdata), e))
            raise  # usually don't need to raise, will kill the script. Re-enable
                   #if need to debug and find where errors are
    print('Time taken : {}\n'.format(datetime.now() - startTime))
    return 0

    # ============== Load refseq and start matching peptides ================ #
    logging.info('Time taken : {}.\n\n'.format(datetime.now() - startTime))
    if not usedb:
        return failed_exps
    else:
        if failed_exps:  # list of failed experiments, is empty if no failures

            for failed in failed_exps:  # failed is a tuple with the datafile and
                                        #then the reason for failure
                usrdatas.remove(failed[0])  # remove all of the failed
                                            #experiments so they won't go in log
                conn = ispec.filedb_connect()
                cursor = conn.cursor()
                cursor.execute("""UPDATE iSPEC_ExperimentRuns
                SET exprun_Grouper_FailedFLAG = ?
                WHERE exprun_EXPRecNo= ? AND
                exprun_EXPRunNo = ? AND
                exprun_EXPSearchNo = ?
                """, 1, exp_setup['EXPRecNo'],
                                exp_setup['EXPRunNo'], exp_setup['EXPSearchNo'])
                conn.commit()
                failedlog.write('{} : failed grouping {},'\
                                ' reason : {}\n'.format(datetime.now(), *failed))
