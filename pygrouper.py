#===============================================================================#
# PyGrouper - Alex Saltzman
import re, csv, os, sys, itertools, time
import logging
import argparse
import threading
from collections import namedtuple, defaultdict
from math import ceil
from statistics import mean
# from sets import Set #python2 only
from configparser import ConfigParser, ParsingError, NoSectionError, NoOptionError
import pandas as pd
from subfuncts import *
try :
    import database_config as db
except ImportError:
    print('Not using databse_config')
    pass # don't use db


program_title = 'PyGrouper v0.1.010'
release_date = '10 December 2015'
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
logfilename = program_title.replace(' ', '_') + '.log'
logging.basicConfig(filename=logfilename, level=logging.DEBUG)
logging.info('{}: Initiating {}'.format(datetime.now(), program_title))

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
    pickle.dump( df, open( name, 'wb' ) )
    print('Pickling...')
    if q:
        print('Exiting prematurely')
        sys.exit(0)

def PeptidomeMatcher(usrdata, ref_dict):
    usrdata['metadatainfo'] = usrdata.apply(lambda x:
                                            genematcher(x['Sequence'],
                                                        x['metadatainfo'],
                                                        ref_dict), axis=1)
    return usrdata


def Grouper(usrfile, usrdata, exp_setup, FilterValues, usedb=False, outdir='', *args):
    #import RefseqInfo

    global program_title
    global release_date
    usrfile = os.path.split(usrfile)[1]
    # file with entry of gene ids to ignore for normalizations

    gid_ignore_file = 'pygrouper_geneignore.txt'
    gid_ignore_list = []

    if os.path.isfile(gid_ignore_file):
        print('Using gene filter file for normalization.')
        gid_ignore_read=[x.strip() for x in open(gid_ignore_file,'r') if
                         not x.strip().startswith('#')]
        #gid_ignore_list = [int(x) for x in gid_ignore_read if x.isdigit()]
        # don't convert to int, GIDs are not ints
        gid_ignore_list = gid_ignore_read

    if exp_setup['EXPQuantSource'].strip() == 'AUC':
        area_col = 'PrecursorArea'  # we always use Precursor Area
        normalize = 10**9
    elif exp_setup['EXPQuantSource'].strip() == 'Intensity':
        area_col = 'Intensity'
        normalize = 10**5    

                
    print('Starting Grouper for exp file {}'.format(usrfile))
    logfilestr = '_'.join(str(x) for x in [exp_setup['EXPRecNo'],
                                           exp_setup['EXPRunNo'],
                                           exp_setup['EXPSearchNo'],
                                           exp_setup['EXPTechRepNo'],
                                           str(exp_setup['EXPLabelType']) +
                                           '.log'])
    logfile = open(os.path.join(outdir,logfilestr), 'w+')  # append to previously opened log file
    logfile.write('{} | Starting {} for file : {}\n'.format(
        time.ctime(),program_title, usrfile))

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
        usrdata['psm_TaxonCount'], usrdata['psm_ProteinList'], \
        usrdata['psm_ProteinCount'] = list(zip(
        *usrdata.apply(lambda x : genelist_extractor(x['metadatainfo'],
                                                     ),
                       axis=1)))


    # ==================== Populate gene info ================================ #
    
    # -------------------- Populate gene, protein, homologene lists----------- #
    #usrdata['_data_tGeneList'], usrdata['_data_tProteinList'],\
    #usrdata['_data_GeneCount'], usrdata['_data_ProteinCount'],\
    #usrdata['_data_tHIDList'], usrdata['_data_HIDCount'],\
    #usrdata['_data_ProteinCapacity'], usrdata['_data_tTaxonIDList'],\
    #usrdata['_data_TaxonCount'] = list(zip(
    #    *usrdata.apply(lambda x:
    #                   meta_extractor(x['metadatainfo']),axis=1)))
    #------------------------------------------------------------------------- #

    usrdata['psm_HID'], usrdata['psm_ProteinGI'] = '', ''
    
    # potentially filled in later,
    #fields will exist in database at least

    logfile.write('{} | Finished matching PSMs to {} '\
                  'refseq\n'.format(time.ctime(),
                                    usrdata.loc[0]['psm_TaxonID']))
    logging.info('Starting Grouper for exp number'\
                 '{}'.format(exp_setup['EXPRecNo']))
        
    nomatches = usrdata[usrdata['psm_GeneCount'] == 0].Sequence  # select all
    #PSMs that didn't get a match to the refseq peptidome
    logfile.write('{} | Total identified PSMs : {}\n'.format(
        time.ctime(), len(usrdata[usrdata['psm_GeneCount'] != 0])))
    logfile.write('{} | Total unidentified PSMs : {}\n'.format(time.ctime(),
                                                               len(nomatches)))
    for missing_seq in nomatches:
        # print('No match for sequence {} in {}').format(missing_seq,ursfile[0])
        # No need to print this to screen
        logging.warning('No match for sequence {} in {}'.format(missing_seq,
                                                                usrfile))
        # Store all of these sequences in the big log file, not per experiment.

    resultout_dir = os.getcwd()  

    logfile.write('{} | Starting grouper.\n'.format(time.ctime()))
    usrdata['psm_PSM_IDG'] = usrdata.apply(lambda x:
                                             IDG_picker(x['IonScore'],
                                                        x['q_value']), axis=1) 
    #usrdata['Sequence'], usrdata['_data_SequenceModi'],
    #usrdata['_data_SequenceModiCount'], usrdata['_data_LabelFLAG'] =
    #list(zip(*usrdata.apply(lambda x : \
        #seq_modi(x['Sequence'], x['Modifications']), axis=1)))
    # We now move this to earlier, before we search against refseq
                                  #(this will deal with the Xs)
    # ============= Gather all genes that all peptides map to =============== #
    glstsplitter = usrdata['psm_GeneList'].str.split(',').apply(pd.Series,
                                                                   1).stack()
    glstsplitter.index = glstsplitter.index.droplevel(-1)  # get rid of
    # multi-index
    glstsplitter.name = 'psm_GeneID'  # give the series a name
    usrdata = usrdata.join(glstsplitter)  # usrdata gains column 'psm_GeneID'
                                          #from glstsplitter Series

    # ========================================================================= #

    logfile.write('{} | Starting peptide ranking.\n'.format(time.ctime()))
    usrdata = usrdata.sort_values(by=['SpectrumFile', 'psm_GeneID', area_col,
                            'Sequence', 'Modifications',
                                      'Charge','psm_PSM_IDG','IonScore', 'PEP',
                            'q_value'], ascending=[0, 1, 0, 1, 1, 1, 1, 0, 1, 1]) 
    usrdata.reset_index(inplace=True)
    usrdata.Modifications.fillna('', inplace=True)  # must do this to compare nans
    usrdata[area_col].fillna(0, inplace=True)  # must do this to compare
                                                       #nans
    grouped = usrdata.drop_duplicates(subset=['SpectrumFile', 'psm_GeneID',
                                              'Sequence', 'Modifications',
                                              'Charge', area_col]).groupby(
                                                  ['SpectrumFile',
                                                   'psm_GeneID'])  # each group
    #belongs to 1 gene, now we can rank on a per-gene basis
    ranks = grouped.cumcount() + 1  # add 1 to start the peptide rank at 1, not 0
    ranks.name = 'psm_PeptideRank'
    usrdata = usrdata.join(ranks)
    usrdata['psm_PeptideRank'] = usrdata['psm_PeptideRank'].fillna(0)  # anyone who
                              # didn't get a rank gets a rank of 0
    logfile.write('{} | Peptide ranking complete.\n'.format(time.ctime()))
    print('{}: Peptide ranking complete for {}.'.format(datetime.now(), usrfile))
    logging.info('{}: Peptide ranking complete for {}.'.format(datetime.now(),
                                                               usrfile))
    # ========================================================================= #
    #quick_save(usrdata, q=False)
    #quick_save(gene_metadata, name='metadata.p', q=False)
    
    usrdata['psm_AUC_useflag'], usrdata['psm_PSM_useflag'] = \
    list(zip(*usrdata.apply(AUC_PSM_flagger, args=(FilterValues,), axis=1)))
    # Flag good quality peptides

        # ======================== Plugin for multiple taxons  ===================== #

    usrdata['gene_taxon_map'] = usrdata.apply(lambda x : gene_to_taxon(
        x['psm_GeneID'], gene_taxon_dict), axis=1)
    
    taxon_ids = set(','.join(x for x in usrdata.psm_TaxonIDList.tolist()
                             if x).split(','))    
    #area_col_new = 'psm_Area_taxonAdj'
    #quick_save(usrdata, q=False)
    #quick_save(gene_taxon_dict, name='gene_taxon_dict.p', q=True)
    #if len(taxon_ids) == 1:  # just 1 taxon id present
    #    usrdata[area_col_new] = usrdata[area_col]
    if len(taxon_ids) > 1:  # more than 1 taxon id
        print('Multiple taxons found, redistributing areas...')
        logfile.write('{} | Multiple taxons found, '\
                      'redistributing areas.\n'.format(time.ctime()))
        usrdata.reset_index(inplace=True)
        #usrdata[area_col_new] = 0
        taxon_totals = dict()
        for taxon in taxon_ids:
            #all_others = [x for x in taxon_ids if x != taxon]
            uniq_taxon = usrdata[
                #(usrdata._data_tTaxonIDList.str.contains(taxon)) &
                #(~usrdata._data_tTaxonIDList.str.contains('|'.join(all_others)))&
                (usrdata['psm_TaxonIDList']==taxon) &
                #(usrdata['psm_PSM_IDG']<9) &  # this is redunant with AUC_UseFLAG
                (~usrdata['psm_GeneID'].isin(gid_ignore_list)) &
                (usrdata['psm_AUC_useflag'] == 1)
                ]
            taxon_totals[taxon] = (uniq_taxon[area_col] / uniq_taxon['psm_GeneCount']).sum()


        tot_unique = sum(taxon_totals.values())  #sum of unique
        # now compute ratio:
        for taxon in taxon_ids:
            taxon_totals[taxon] = taxon_totals[taxon] / tot_unique
            print(taxon, ' ratio : ', taxon_totals[taxon])
            logfile.write('{} ratio : {}\n'.format(taxon, taxon_totals[taxon]))

        ### We don't want to do this... ###
        #all_combos = [x for i in range(2, len(taxon_ids)+1) for x in
        #              itertools.combinations(taxon_ids, i)] # list of tuples
        #patterns = regex_pattern_all(all_combos)
        #for taxons, pattern in zip(all_combos, patterns):
        #    for taxon in taxons:
        #        ratio = taxon_totals[taxon]
        #        usrdata.ix[(usrdata.psm_tTaxonIDList.str.contains(pattern)) &
        #                   (usrdata.gene_taxon_map == taxon),
        #                   area_col_new] = usrdata[area_col] * ratio
        #
        #usrdata.ix[usrdata.psm_TaxonCount==1, area_col_new] = usrdata[area_col]
        #area_col = area_col_new  # use new area col as the area column now
        print()
        #sys.exit(0)    
    # ========================================================================= #
            

    
    pd.options.mode.chained_assignment = None  # default='warn'

    # Make the name for the peptide data table :
    usrdata_out = '_'.join(str(x) for x in [exp_setup['EXPRecNo'],
                                            exp_setup['EXPRunNo'],
                                            exp_setup['EXPSearchNo'],
                                            exp_setup['EXPTechRepNo'],
                                            exp_setup['EXPLabelType'],
                                            'data.tab'])
    # none/SILAC loop
    labeltypes = ['nolabel', 'heavy']  # right now only using nolabel, but in
                                       # future this can grow
    gpgcount, genecount, ibaqtot, = 0, 0, 0                                   
    for label in range(1):  # increase the range to go through more label types
        logfile.write('{} | Starting gene assignment for label {}.\n'.format(
            time.ctime(), labeltypes[label]))
        # ==========Select only peptides flagged  with good quality=========== #
        temp_df = usrdata[(usrdata['psm_LabelFLAG'] == label) &
                          (usrdata['psm_AUC_useflag'] == 1) &
                          (usrdata['psm_GeneCount'] > 0)]  # should keep WL's
                                            #peptides too (if they have a mass)
        gene_cols = ['gene_EXPRecNo', 'gene_EXPRunNo', 'gene_EXPSearchNo',
                     'gene_EXPLabelFLAG', 'gene_AddedBy',
                     'gene_CreationTS', 'gene_ModificationTS', 'gene_GeneID',
                     'gene_IDSet', 'gene_IDGroup','gene_IDGroup_u2g',
                     'gene_GPGroup', 'gene_GPGroups_All', 'gene_PSMs',
                     'gene_PSMs_u2g', 'gene_PeptidePrint', 'gene_PeptideCount',
                     'gene_PeptideCount_u2g', 'gene_PeptideCount_S',
                     'gene_PeptideCount_S_u2g','gene_GeneArea_gpcAdj',
                     'gene_GeneArea_gpcAdj_u2g', 'gene_GeneArea_gpcAdj_u2g_all',
                     'gene_GeneArea_gpcAdj_max', 'gene_GeneArea_dstrAdj',
                     'gene_GeneCapacity', 'gene_iBAQ']  # cols of
                                                                  #interest

        # ==================================================================== #
        #print(len(temp_df))  # for debugging
        if len(temp_df) > 0:  # only do if we actually have peptides selected
            genedata_out = '_'.join(str(x) for x in [exp_setup['EXPRecNo'],
                                                     exp_setup['EXPRunNo'],
                                                     exp_setup['EXPSearchNo'],
                                                     exp_setup['EXPTechRepNo'],
                                                     exp_setup['EXPLabelType'],
                                                     labeltypes[label],
                                                     'e2g.tab'])
            print('{}: Populating gene table for {}.'.format(datetime.now(),
                                                             usrfile))
            # genes_df['_e2g_GeneID'] = Set(temp_df['_data_GeneID']) #only in 2.7
            genes_df = pd.DataFrame({'gene_GeneID':
                                     list(set(temp_df['psm_GeneID']))})

            genes_df['gene_GeneCapacity'] = genes_df.apply(lambda x:
                                                           capacity_grabber(
                                                               x['gene_GeneID'],
                                                               gene_metadata),
                                                           axis=1)
            genes_df['gene_EXPLabelFLAG'] = label
            #for k, v in exp_setup.items():
            #    if k.startswith('EXP'):
            #        genekey = [key for key in exp_setup if
            #                   key.startswith('gene_') and k[3:] in key]
            #        print(genekey)
            #        if len(genekey)==1:
            #            genekey = genekey[0]
            #            genes_df[genekey] = v  # creates genes_df
            #columns and populates them with  exp_setup values.
            #Quick, easy, labeling

            genes_df['gene_PeptideSet'], genes_df['gene_PeptidePrint'], \
            genes_df['gene_PeptideCount'], genes_df['gene_PeptideCount_u2g'],\
            genes_df['gene_PeptideCount_S'],\
            genes_df['gene_PeptideCount_S_u2g'] =\
            list(zip(*genes_df.apply(pept_print, args=(temp_df,), axis=1)))

            genes_df['gene_PSMs'], genes_df['gene_PSMs_u2g'],\
            genes_df['gene_PSMs_S'],genes_df['gene_PSMs_S_u2g'] = \
            list(zip(*genes_df.apply(lambda x:
                                     e2g_PSM_helper(x['gene_GeneID'],
                                                    temp_df,
                                                    exp_setup['EXPTechRepNo']),
                                     axis=1)))

            print('{}: Calculating peak areas for {}.'.format(
                datetime.now(), usrfile))
            logging.info('{}: Calculating peak areas for {}.'.format(
                datetime.now(), usrfile))
            logfile.write('{} | Calculating peak areas.\n'.format(time.ctime()))

            genes_df['gene_GeneArea_gpcAdj_max'],genes_df['gene_GeneArea_gpcAdj'],\
            genes_df['gene_GeneArea_gpcAdj_u2g'],\
            genes_df['gene_GeneArea_gpcAdj_u2g_all']  = \
            list(zip(*genes_df.apply(area_calculator,
                                     args=(temp_df,
                                           exp_setup['EXPTechRepNo'],
                                           area_col, normalize),
                                     axis=1)))
            # pandas may give a warning from this though it is fine
            print('{}: Calculating distributed area ratio for {}.'.format(
                datetime.now(), usrfile))
            logging.info('{}: Calculating distributed area ratio for {}.'.format(
                datetime.now(), usrfile))
            logfile.write('{} | Calculating distributed area ratio.\n'.format(
                time.ctime()))
            #quick_save(genes_df,name='genesdf_snapshot.p', path=None, q=False)
            #quick_save(temp_df,name='tempdf_snapshot.p', path=None, q=True)
            temp_df['psm_PrecursorArea_dstrAdj'] = temp_df.apply(
                AUC_distributor,args=(genes_df,
                                      area_col,),
                                      axis=1)
            
            print('{}: Assigning gene sets and groups for {}.'.format(
                datetime.now(), usrfile))
            logging.info('{}: Assigning gene sets and groups for {}.'.format(
                datetime.now(), usrfile))
            logfile.write('{} | Assigning gene sets and groups.\n'.format(
                time.ctime()))

            #quick_save(genes_df,name='genes_df_snapshot.p', path=None, q=False)
            #quick_save(temp_df,name='temp_df_snapshot.p', path=None, q=True)
            genes_df['gene_IDSet'], genes_df['gene_IDGroup'],\
            genes_df['gene_IDGroup_u2g'] = list(
                zip(*genes_df.apply(gene_setter, args=(
                    genes_df, temp_df,), axis=1)))

            genes_df['gene_GeneArea_dstrAdj'] = genes_df.apply(gene_AUC_sum,
                                                                  args=(temp_df,
                                                                  normalize,), 
                                                                  axis=1)
            genes_df['gene_GeneCapacity'] =\
                                        genes_df.gene_GeneCapacity.astype('float')
            # print(genes_df._e2g_nGPArea_Sum_dstrAdj.dtype,1
            #genes_df._e2g_GeneCapacity.dtype)  # for debugging
            genes_df['gene_iBAQ'] = \
                genes_df.gene_GeneArea_dstrAdj / genes_df.gene_GeneCapacity
            genes_df['gene_GPGroup'] = ''
            genes_df.sort_values(by=['gene_PSMs'], ascending=False, inplace=True)
            genes_df.index = list(range(0, len(genes_df)))
            last = 0
            #quick_save(genes_df,name='genes_df_snapshot.p', path=None, q=False)
            #quick_save(temp_df,name='temp_df_snapshot.p', path=None, q=True)
            for i in range(len(genes_df)):  # The logic behind it makes
                #sense,
                #but the implementation is weird.
            # print(last)  # for debugging
                if genes_df.loc[i]['gene_IDSet'] != 3:
                    genes_df.loc[i,'gene_GPGroup'], lessthan = \
                    GPG_helper(genes_df.at[i,'gene_IDSet'],
                               genes_df.at[i,'gene_PeptideSet'], \
                               genes_df, last)
                if type(genes_df.loc[i]['gene_GPGroup']) is int and not lessthan:
                    last = genes_df.loc[i, 'gene_GPGroup']

            genes_df['gene_GPGroups_All'] = genes_df.apply(GPG_all_helper,
                                                           args=(genes_df,),
                                                           axis=1)
            genes_df['gene_GPGroup'].replace(to_replace='', value=float('NaN'),
                                         inplace=True)  # can't sort int and
                                         #strings, convert all strings to NaN
            genes_df.sort_values(by=['gene_GPGroup'], ascending=True, inplace=True)
            genes_df.index = list(range(0, len(genes_df)))  # reset the index
            gpgcount += genes_df.gene_GPGroup.max()  # do this before filling na
            genes_df['gene_GPGroup'].fillna('', inplace=True)  # convert all NaN
                                                           #back to empty string
            # =============================================================#
            genes_df['gene_EXPRecNo'] = exp_setup['EXPRecNo']
            genes_df['gene_EXPRunNo'] = exp_setup['EXPRunNo']
            genes_df['gene_EXPSearchNo'] = exp_setup['EXPSearchNo']
            genes_df['gene_EXPTechRepNo'] = exp_setup['EXPTechRepNo']
            genes_df['gene_AddedBy'] = usrdata.loc[1]['psm_AddedBy']
            genes_df['gene_CreationTS'] = datetime.now().ctime()
            genes_df['gene_ModificationTS'] = datetime.now().ctime()
            #quick_save(genes_df,name='genes_df_snapshot.p', path=None, q=False)
            #quick_save(gene_metadata,name='genemetadata_snapshot.p', path=None, q=False)
            #quick_save(temp_df,name='tempdf_snapshot.p', path=None, q=False)
            renamed_genecols = [exp_setup.get(genecol, genecol) for genecol in gene_cols]
            torename = {k:v for k,v in exp_setup.items() if k.startswith('gene_')}
            # get a log of gene info
            genecount += len(genes_df)
            ibaqtot += genes_df[~genes_df.gene_GeneID.isin(
                gid_ignore_list)].gene_iBAQ.sum()
            if usedb:
                #db.add_exp2gene(genes_df)
                pass

            genes_df.rename(columns=torename, inplace=True)
            genes_df.to_csv(os.path.join(outdir,genedata_out), columns=renamed_genecols,
                            index=False,encoding='utf-8', sep='\t')
            logfile.write('{} | Export of genetable for labeltype {}'\
                          'completed.\n'.format(
                              time.ctime(), 
                              labeltypes[label]))


    # ----------------End of none/silac loop--------------------------------- #
    usrdata.drop('metadatainfo', axis=1, inplace=True)  # Don't need this
                                      # column anymore.

    
    usrdata = pd.merge(usrdata, temp_df, how='outer')
    usrdata['psm_EXPRecNo'], usrdata['psm_EXPRunNo'],\
    usrdata['psm_EXPSearchNo'],\
    usrdata['psm_EXPTechRepNo'] = exp_setup['EXPRecNo'],\
    exp_setup['EXPRunNo'], exp_setup['EXPSearchNo'],\
    exp_setup['EXPTechRepNo']

    usrdata['psm_CreationTS'] = datetime.now().ctime()
    usrdata['psm_ModificationTS'] = datetime.now().ctime()
    usrdata['psm_HID_list'] = ''  # will be populated later
    usrdata['psm_HID_count'] = ''  # will be populated later
    data_cols = ['psm_EXPRecNo', 'psm_EXPRunNo', 'psm_EXPSearchNo',
                 'psm_EXPTechRepNo', 'Sequence',
                 'PSMAmbiguity', 'Modifications', 'ActivationType',
                 'DeltaScore', 'DeltaCn', 'Rank', 'SearchEngineRank',
                 'PrecursorArea', 'QuanResultID', 'q_value', 'PEP',
                 'Decoy Peptides Matched', 'Exp Value', 'Homology Threshold', 
                 'Identity High', 'Identity Middle', 'IonScore',
                 'Peptides Matched', 'MissedCleavages',
                 'IsolationInterference', 'IonInjectTime',
                 'Intensity', 'Charge', 'mzDa', 'MHDa',
                 'DeltaMassDa', 'DeltaMassPPM', 'RTmin',
                 'FirstScan', 'LastScan', 'MSOrder', 'MatchedIons', 
                 'TotalIons', 'SpectrumFile', 'Annotation', 'psm_AddedBy',
                 'psm_CreationTS','psm_ModificationTS', 'psm_GeneID',
                 'psm_GeneList', 'psm_GeneCount', 'psm_ProteinGI',
                 'psm_ProteinList','psm_ProteinCount', 
                 'psm_HID', 'psm_HID_list', 'psm_HID_count',
                 'psm_TaxonID', 'psm_TaxonIDList', 'psm_TaxonCount',
                 'psm_PSM_IDG', 'psm_SequenceModi', 
                 'psm_SequenceModiCount', 'psm_LabelFLAG', 
                 'psm_PeptideRank', 'psm_AUC_useflag','psm_PSM_useflag',
                 'psm_PrecursorArea_dstrAdj']
    #usrdata.to_csv(usrdata_out, columns=usrdata.columns,
                                #encoding='utf-8', sep='\t')
    #print(usrdata.columns.values)  # for debugging
    if not all(x in usrdata.columns.values for x in data_cols):
        print('Potential error, not all columns filled.')
        print([x for x in data_cols if x not in usrdata.columns.values])                                
        data_cols = [x for x in data_cols if x in usrdata.columns.values]
        # will still export successfully

    # update database
    if usedb:
        session = db.make_session()
        exprecord = session.query(db.ExperimentRun).\
                    filter_by(record_no=exp_setup['EXPRecNo']).\
                    filter_by(run_no=exp_setup['EXPRunNo'],
                              search_no=exp_setup['EXPSearchNo'],
                              tech_repeat=exp_setup['EXPTechRepNo']).one()
        exprecord.GPGroup_count = int(gpgcount)
        exprecord.gene_count = int(genecount)
        exprecord.iBAQ_total = int(ibaqtot)
        exprecord.PSM_count = len(usrdata[usrdata.psm_PSM_useflag==1])
        exprecord.grouped = True
        exprecord.failed = False
        exprecord.group_date = datetime.now()

        session.add(exprecord)
        session.commit()
        session.close()

    renamed_datacols = [exp_setup.get(datacol, datacol) if datacol.startswith('psm_') else datacol
                        for datacol in data_cols]
    usrdata.rename(columns={k:v for k,v in exp_setup.items() if k.startswith('psm_')}, inplace=True)
    usrdata.to_csv(os.path.join(outdir,usrdata_out), columns=renamed_datacols,
                   index=False, encoding='utf-8',sep='\t')
    
    logfile.write('{} | Export of datatable completed.\n'.format(time.ctime()))
    logfile.write('Successful grouping of file completed.')
    logfile.close()
    

                                
    print('Successful grouping of {} completed.\n' \
          .format('_'.join(
              [str(exp_setup['EXPRecNo'])+'_'+str(exp_setup['EXPRunNo']
              )])))

def main(usrfiles=[], exp_setups=[], automated=False, usepeptidome=True, setup=False,fullpeptread=False,
         usedb=False, inputdir='', outputdir=''):
    # ===================Configuration Setup / Loading==========================#
    global program_title
    global release_date
    refs = {}
    if not os.path.isfile(os.path.join(BASE_DIR,'py_config.ini')):
        input("No config file detected. Don't worry, we'll make one now\n"\
              "Press [Enter] to continue")
        pysetup()
    elif setup:
        pysetup()
    parser = ConfigParser(comment_prefixes=(';')) # allow number sign to be read in configfile
    parser.optionxform = str  # preserve case
    parser.read(os.path.join(BASE_DIR,'py_config.ini'))
    pept_breakups = {}
    breakup_size = 4
    for taxon, location in parser.items('refseq locations'):
        refs[taxon] = {'loc': location,
                       'size': parser.getfloat('refseq file sizes',
                        taxon)}  # access and store preconfigured reference info
        if fullpeptread:
            pept_breakups[taxon] = [(0, int(refs[taxon]['size']))]
        else:
            pept_breakups[taxon] = [(0, int(ceil(refs[taxon]['size'] / \
                                                 breakup_size)))] * breakup_size

    # ====================Configuration Setup / Loading======================= #

    if imagetitle:
        fancyprint(program_title, 12)  # ascii art
        #  fancyprint('Malovannaya lab',10)  #
    elif not imagetitle:
        print(program_title)  # not as exciting as ascii art

    print('\nrelease date: {}'.format(release_date))
    print('Python version ' + sys.version)
    print('Pandas version: ' + pd.__version__)

    if not automated: usr_name = input('Enter your name : ')

    FilterValues = {'Filter_IS': 7, 'Filter_qV': 0.05, 'Filter_PEP': 'all',
                    'Filter_IDG': 'all', 'Filter_Z_min': 2,'Filter_Z_max': 4,
                    'Filter_Modi': 3}  # defaults
    Filter_Stamp ='is{}_qv{}_pep{}_idg{}_z{}to{}mo{}'.format(
        FilterValues['Filter_IS'],
        FilterValues['Filter_qV'] * 100,
        FilterValues['Filter_PEP'],
        FilterValues['Filter_IDG'], 
        FilterValues['Filter_Z_min'],
        FilterValues['Filter_Z_max'], 
        FilterValues['Filter_Modi'])

    print('\nFilter values set to : {}'.format(Filter_Stamp))
    if not automated and not usedb:
        try:
            explog = open('PyGrouper_grouped_exps.log', 'r+U')
        except IOError:
            explog = open('PyGrouper_grouped_exps.log', 'a+')
        grouped_exps = [value for value in re.findall(r'(\d+\S*\.txt)',\
                                                    ' '.join(explog.readlines()))]
    if not automated and usedb:
        grouped_query = db.get_grouped_exps()
        grouped_exps = [(group.record_no, group.run_no, group.search_no, group.tech_repeat)
                        for group in grouped_query]
    if not automated:
        try:
            input('Press enter to continue, or Ctrl+C to modify the filter'\
                  'values.\n')
        except KeyboardInterrupt:
            print('\nFilter_PEP and Filter_IDG can be set to "all",'\
                  ' but all other values must be numbers only.')
            dict_modifier(FilterValues, {'Filter_PEP': ['all'],
                                         'Filter_IDG': ['all'], 'Filter_qV': []})

        logging.info('Filter stamp : {}'.format(Filter_Stamp))

        while True:
            try:
                exp_setup = {'taxonID': 9606, 'EXPTechRepNo': 1,
                             'EXPQuantSource': 'AUC', 'EXPRunNo': 1,
                             'EXPSearchNo': 1, 'EXPLabelType': 'none'}
                usrfile_input = input('Enter a file to group or press'\
                                      ' Ctrl+C if done : ')
                if usrfile_input == 'forcequit':
                    logging.info('forcequit without selecting any files')
                    logging.shutdown()
                    sys.exit(0)
                if os.path.isfile(usrfile_input):  # check to see if file exists
                    proceed = True
                    if usrfile_input.strip() in grouped_exps and not usedb:  # strip any
                        #whitespace to match correctly
                        while True:
                            proceed = input(
                                'Records show that {} has been grouped before.'\
                                'Are you sure you want to regroup (Y/n)? '\
                                .format(usrfile_input))
                            if 'Y' in proceed:
                                break
                            elif 'n' in proceed.lower():
                                proceed = False
                                break
                    if proceed:
                        try:
                            exp_setup['EXPRecNo'] = int(re.search('(\d{3,})',
                                                    usrfile_input).group())
                            # find first number of at least 3 digits
                        except AttributeError:
                            exprecno = input("Couldn't locate experimental"\
                                             " record automatically,"\
                                             " please enter it now. ")
                            exp_setup['EXPRecNo'] = int(exprecno)

                        print('Experimental setup is : {}'.format(exp_setup))
                        try:
                            input('Press enter to accept values and continue,'\
                                  'or Ctrl+C to modify the experimental'\
                                  'values.\n')
                        except KeyboardInterrupt:
                            print(
                                '\nEXPQuantSource can be set to AUC or '\
                                'Intensity, but all other values must be '\
                                'numbers only.')
                            dict_modifier(exp_setup, {'EXPQuantSource':
                                                      ['AUC', 'Intensity']}, 
                                          exception_float=False)
                            
                        if usedb:
                            #print(grouped_exps)
                            #print((exp_setup['EXPRecNo'], exp_setup['EXPRunNo'],
                            #    exp_setup['EXPSearchNo'], exp_setup['EXPTechRepNo']))
                            exp_setup['add_to_db'] = True
                            if (exp_setup['EXPRecNo'], exp_setup['EXPRunNo'],
                                exp_setup['EXPSearchNo'], exp_setup['EXPTechRepNo']) in grouped_exps:
                                while True:
                                    proceed = input(
                                        'Records show that {} has been grouped before.'\
                                        'Are you sure you want to regroup (Y/n)? '\
                                        .format(usrfile_input))
                                    if 'Y' in proceed:
                                        exp_setup['add_to_db'] = False
                                        break
                                    elif 'n' in proceed.lower():
                                        proceed = False
                                        break
                        if proceed:            
                            exp_setups.append(exp_setup)
                            usrfiles.append(usrfile_input)
                                                        
                else:
                    print('File {} not found.'.format(usrfile_input))
                    # morefile = raw_input('Do you have more files to group? ')
            except KeyboardInterrupt:
                if len(usrfiles) > 0:
                    print()
                    if usedb:
                        newexps = defaultdict(list)
                        for exp in exp_setups:
                            if exp['add_to_db']:
                                newexps[exp['EXPRecNo']].append(
                                    {'run_no':exp['EXPRunNo'],
                                     'search_no':exp['EXPSearchNo'],
                                     'taxon':exp['taxonID'],
                                     'addedby':usr_name,
                                     'creation_ts':datetime.now(),
                                     'techrep':exp['EXPTechRepNo'],
                                     'label':exp['EXPLabelType'],
                                     'quant':exp['EXPQuantSource'],
                                     })
                            
                        db.add_experiments(newexps)
                    break
                else:
                    print('No files selected!')

    startTime = datetime.now()
    if FilterValues['Filter_PEP'] == 'all' or FilterValues['Filter_IDG'] == 'all':
        FilterValues['Filter_PEP'] = float('inf')

    print('\nStart at {}'.format(startTime))
    logging.info('Start at {}'.format(startTime))
    if fullpeptread: print('Running with fullpeptread option')
    usrdatas = []
    column_aliases = dict(parser.items('column names'))
    for key in column_aliases:  # find the right column name
        if key.startswith('psm_') or key.startswith('gene_'):
            column_aliases[key] = list(filter(None,
                                              (x.strip() for
                                               x in column_aliases[key].splitlines())))[-1]
            
        else:
            column_aliases[key] = list(filter(None,
                                          (x.strip() for
                                           x in column_aliases[key].splitlines())))
    for usrfile in usrfiles:  # load all data
        usrdatas.append(pd.read_csv(os.path.join(inputdir,usrfile), sep='\t'))
    for usrdata, exp_setup in zip(usrdatas, exp_setups):
        standard_names = column_identifier(usrdata, column_aliases)
        for name in standard_names:
            exp_setup[name] = standard_names[name]
        for alias in column_aliases:
            if alias.startswith('psm_') or alias.startswith('gene_'):
                exp_setup[alias] = column_aliases[alias]

        #for key in exp_setup: print(key,'   :   ', exp_setup[key])
        usrdata.rename(columns={v: k
                                for k,v in standard_names.items()},
                       inplace=True)
        
        
        #usrdata.rename(columns={'Annotated Sequence': 'Sequence', 'Area':
        #                        'Precursor Area','Percolator q-Value':
        #                        'q-Value','Percolator PEP': 'PEP',
        #                        'Ions Score': 'IonScore', 'DeltaM [ppm]':
        #                        'Delta Mass [PPM]',
        #                        'Deltam/z [Da]' : 'Delta Mass [Da]',
        #                        }, inplace=True)  # PD 2.0 column name changes
        if not automated:
            usrdata['psm_AddedBy'] = usr_name

        elif automated:
            usrdata['psm_AddedBy'] = exp_setup['AddedBy']

        usrdata['psm_TaxonID'] = exp_setup['taxonID']
        usrdata['psm_GeneList'], usrdata['psm_ProteinList'],\
        usrdata['psm_GeneCount'],usrdata['psm_ProteinCount'],\
        usrdata['psm_HomologeneID'], usrdata['psm_ProteinCapacity'], \
        usrdata['metadatainfo'] = '', '', 0, 0, '', '', ''

    # ============== Load refseq and start matching peptides ================ #
    print('{}: Loading refseq database.'.format(datetime.now()))
    #RefseqInfo = namedtuple('RefseqInfo',
    #                    'taxonid, geneid, homologeneid,proteingi,genefraglen')

    for organism in refs.keys():
        if any(any(x['psm_TaxonID'].isin([int(organism)])) for x in\
               usrdatas):  # check if we even need to read the
                           #peptidome for that organism
            ref_reader = csv_reader(refs[organism]['loc'])
            print('Using peptidome {} '.format(refs[organism]))
            #print('Breakups : {}'.format(pept_breakups[organism]))
            #sys.exit(0)
            for breakup in pept_breakups[organism]:  # Read refseq in chunks,
                                                     #uses less memory
                prot = defaultdict(list)
                for k in range(*breakup):
                    try:
                        row = next(ref_reader)
                        fragments, fraglen = protease(row.fasta, minlen=7,
                                                      cutsites=['K', 'R'],
                                                      exceptions=['P'])
                        # should stop else clause here (??)
                        # fragments = protease((linesplit[4].strip('\n'),minlen =
                                                  #7, cutsites=['k','r'])
                        for fragment in fragments:
                            prot[fragment].append(
                                RefseqInfo._make([row.taxonid, row.geneid,
                                                  row.homologeneid,row.proteingi, 
                                                  fraglen]))

                    except StopIteration:  # breakups won't be exact since they
                        #are rounded up to ensure full file coverage ##should
                        #move this up and use else clause
                        break

                for usrdata, usrfile in zip(usrdatas, usrfiles):
                    #print(usrdata.loc[0]['_data_TaxonID'])
                    if usrdata.loc[0]['psm_TaxonID'] == int(organism):
                        # check to see if the inputdata 
                        #taxonID matches with the proteome
                        usrdata['Sequence'], usrdata['psm_SequenceModi'],\
                        usrdata['psm_SequenceModiCount'],\
                        usrdata['psm_LabelFLAG'] = \
                        list(zip(*usrdata.apply(lambda x :
                                                seq_modi(x['Sequence'],
                                                         x['Modifications']),
                                                axis=1)))
                        # print 'Matching for {}'.format(usrfile)
                        PeptidomeMatcher(usrdata, prot)  # call matcher
                del prot  # free up memory

    print('{}: Finished matching peptides to genes.'.format(datetime.now()))
    logging.info('{}: Finished matching peptides to'\
                 'genes.'.format(datetime.now()))
    failed_exps = []
    for usrfile, udata, esetup in zip(usrfiles, usrdatas, exp_setups):
        # print usrfile, esetup
        try:
            Grouper(usrfile, udata, esetup, FilterValues, usedb=usedb, outdir=outputdir)
        except Exception as e:  # catch and store all exceptions, won't crash
                                # the whole program at least
            failed_exps.append((usrfile, e))
            print('Failure for file of experiment {}.\n'\
                  'The reason is : {}'.format(esetup['EXPRecNo'], e))
            logging.warn('Failure for file of experiment {}.\n'\
                         'The reason is : {}'.format(esetup['EXPRecNo'], e))
            raise  # usually don't need to raise, will kill the script. Re-enable
                   #if need to debug and find where errors are
    print('Time taken : {}\n'.format(datetime.now() - startTime))
    logging.info('Time taken : {}.\n\n'.format(datetime.now() - startTime))
    if automated and not usedb:
        return failed_exps
    else:
        if failed_exps:  # list of failed experiments, is empty if no failures

            for failed in failed_exps:  # failed is a tuple with the datafile and
                                        #then the reason for failure
                usrfiles.remove(failed[0])  # remove all of the failed
                                            #experiments so they won't go in log
                if usedb:
                    session = db.make_session()
                    exprecord = session.query(db.ExperimentRun).\
                    filter_by(record_no=exp_setup['EXPRecNo']).\
                    filter_by(run_no=exp_setup['EXPRunNo']).one()
                    exprecord.failed = True
                    session.add(exprecord)
                    session.commit()
                    session.close()
                failedlog.write('{} : failed grouping {},'\
                                ' reason : {}\n'.format(datetime.now(), *failed))
        if usrfile and not usedb:
            for usrfile in usrfiles:
                explog.write('{} : grouped experiment file {}'\
                             '\n'.format(datetime.now(), usrfile))


def expchecker(options):
    '''Checks for new experimental data based on 
    '''
    global program_title
    # exp_setup = {'taxonID' : 9606, 'EXPTechRepNo' : 1, 'EXPQuantSource' :
    # 'AUC', 'EXPRunNo' : 1, 'EXPSearchNo': 9}
    print('{} : Checking for new experiments...'.format(time.ctime()))
    try:
        explog = open('PyGrouper_grouped_exps.log', 'r+U')
    except IOError:  # make a new file if doesn't exist
        explog = open('PyGrouper_grouped_exps.log', 'a+')
    try:
        failedlog = open('PyGrouper_failed_exps.log', 'r+U')
    except IOError:  # make a new file if doesn't exist
        failedlog = open('PyGrouper_failed_exps.log', 'a+')
    cwd = os.getcwd()
    # grouped_exps = [int(value) for value in re.findall(r'(\d+)','
    #'.join(explog.readlines()))]
    grouped_exps = [value for value in re.findall(r'(\d+\S*\.txt)', ' '.join(
        explog.readlines()))]  # starts with at least 1 digit, has unspecified amount of non-whitespace, and ends in .tab
    failed_exps = [value for value in re.findall(r'(\d+\S*\.txt)', \
                                                 ' '.join(failedlog.readlines()))]
    # failed_exps = [int(value) for value in re.findall(r'(\d+)','
    #'.join(failedlog.readlines()))]
    validfiles = [f for f in os.listdir(os.getcwd()) if '.txt' in f]
    ispec_export = '4PyGrouper_ExpRunDump.xlsx'
    experiment_file = pd.read_excel(ispec_export)
    # gen = experiment_file.iterrows()
    exp_setups, usrfiles = [], []
    required_fields = ['_exprun_EXPRecNo','_exprun_EXPRunNo']
    experiment_file.fillna(0, inplace=True)  # Fill NAs with 0 so we can convert
    #everything to int
    for field in required_fields:
        experiment_file[field] = experiment_file[field].astype('int')  # get rid
        #of decimals
    for row in experiment_file.iterrows():  # row is a tuple, [0] is index, [1]
        #is actual row
        if not any(row[1][required_fields].isin([0])):  # ensure all required
            #info is present
            record = str(row[1]._exprun_EXPRecNo) + '_' + \
                     str(row[1]._exprun_EXPRunNo)
            try:
                if all(record not in exp for exp in grouped_exps) and all(
                            record not in fexp for fexp in failed_exps):
                    # then we group
                    exp_setup = {'EXPRecNo': int(row[1]._exprun_EXPRecNo),
                                 'EXPRunNo': int(row[1]._exprun_EXPRunNo),
                                 'EXPSearchNo': int(row[1]._exprun_EXPSearchNo),
                                 'taxonID': int(row[1]._exprun_TaxonID),
                                 'EXPQuantSource':
                                 str(row[1].exprun_Search_QuantSource),
                                 'AddedBy': str(row[1]._exprun_AddedBy),
                                 'EXPTechRepNo': int(row[1].exprun_nTechRepeats),
                                 'EXPLabelType': str(row[1].exprun_LabelType)}
                    passed = True  # flag for if we should group an experiment
                else:
                    passed = False
            except ValueError:
                # print('Experimental information not complete for experiment
                #number')
                passed = False
                pass
        
            if passed:
                usrfilelst = [x for x in validfiles if record in x]
                # match exp record and run numbers
                if len(usrfilelst) == 1:  # ensure we just have 1 file
                    usrfile = usrfilelst[0]
                elif len(usrfilelst) > 1:
                    usrfilelst = [x for x in usrfilelst if 'all' in x]
                    # find the file with 'all' in it,
                    # this is based on how we name files for MudPIT
                    if len(usrfilelst) == 1:
                        usrfile = usrfilelst[0]
                    else:  # can't find the correct file if we get to here
                        print('Warning, more than one match for EXPRecNo {}'\
                              ' detected, skipping...'.format(
                                  exp_setup['EXPRecNo']))
                        passed = False
                elif len(usrfilelst) == 0:  # no file found
                # print 'Warning, no file found for EXPRecNo
                    #{}.'.format(exp_setup['EXPRecNo'])
                    passed = False
            if passed and len(exp_setups) <= 25:  # cap on how may files to group
                                                  #at once, untested
                exp_setups.append(exp_setup)
                usrfiles.append(usrfile)
                passed = False  # reset flag

    if len(usrfiles) > 0:  # then we group
        for usrfile, exp_setup in zip(usrfiles, exp_setups):
            print('Found experiment {} from datafile'\
                  ' {}.'.format(exp_setup['EXPRecNo'], usrfile))
            logfilestr = '_'.join(
                str(x) for x in [exp_setup['EXPRecNo'], exp_setup['EXPRunNo'],
                                 exp_setup['EXPSearchNo'],
                                 exp_setup['EXPTechRepNo'],
                                 exp_setup['EXPLabelType'] + '.log'])
            logfile = open(logfilestr, 'w')
            logfile.write('{} | starting {} for {} with taxonID: {}\n'.format(time.ctime(), program_title, 
                                                                              exp_setup['EXPRecNo'],      
                                                                              exp_setup['taxonID']))
            logfile.close()

        try:
            logging.info('{}: Initiating {}'.format(datetime.now(), program_title))
            failed_exps = main(usrfiles, exp_setups, **options)
            if failed_exps:  # list of failed experiments, is empty if no failures
                for failed in failed_exps:
                    # failed is a tuple with the datafile and then the reason for
                    #failure
                    usrfiles.remove(failed[0])
                    # remove all of the failed experiments so they won't go in log
                    failedlog.write('{} : failed grouping {}, reason :'\
                                    '{}\n'.format(datetime.now(), *failed))
            if usrfile:
                for usrfile in usrfiles:
                    explog.write('{} : grouped experiment file {}'\
                                 '\n'.format(datetime.now(), usrfile))
        except Exception as e:  # catch all exceptions
            logging.exception('Fatal error with Grouper function : {}'.format(e))
            raise  # crash and burn
    failedlog.close()
    explog.close()
    #print('{} : Sleeping... Press Enter to wake or Ctrl+C to
    #end.'.format(time.ctime()))


def schedule(INTERVAL, options):
    expchecker(options)
    global thread
    thread = threading.Timer(INTERVAL, schedule, [INTERVAL, options])
    print('{} : Sleeping... Press Enter to wake or [exit] to'\
          ' end.'.format(time.ctime()))
    thread.start()


if __name__ == '__main__':

    parser = ConfigParser(comment_prefixes=(';')) # allow number sign to be read in configfile
    parser.optionxform = str  # preserve case
    parser.read('py_config.ini')

    INPUT_DIR = parser['directories']['inputdir']  # where to look
                                                   # for files to group
    OUTPUT_DIR = parser['directories']['outputdir']

    parser = argparse.ArgumentParser()
    parser.add_argument('-fpr', '--fullpeptread',
                        action='store_true',
                        help='Default FALSE. Load the peptidome all at once,'
                        'minor speed improvement at the expense of memory. '\
                        'Use at your own risk.')
    parser.add_argument('-gs', '--genesets',
                        action='store_true',
                        help='Optional inclusion of peptide sets for each gene.'\
                        'Legacy and not useful.')
    parser.add_argument('-s', '--setup',
                        action='store_true',
                        help='Run setup wizard for PyGrouper.\n'\
                        'Not necessary unless adding a new refseq or need to'\
                        'change refseq location.')
    parser.add_argument('-a', '--automated',
                        action='store_true', help='Automated run of PyGrouper.'\
                        'Note, requires experiment dump from iSPEC to be set up'\
                        'seperately.\nIf you are not sure if this is correctly'\
                        'set up, automation will probably not work.')
    parser.add_argument('-nd', '--nodatabase', action='store_true', 
                        help='Do not use database to store'\
                        'experiment info.')
    args = parser.parse_args()

    options = {}
    options['setup'] = args.setup
    options['fullpeptread'] = args.fullpeptread
    options['usedb'] = not args.nodatabase
    options['inputdir'] = INPUT_DIR
    options['outputdir'] = OUTPUT_DIR
    if args.automated:
        print('--automated is not obsolete. Please run auto_grouper.py instead')
        sys.exit(0)
        options['automated'] = True
        while True:
            
            INTERVAL = 60 * 60
            schedule(INTERVAL, options)
            #print('{} : Sleeping... Press Enter to wake or [exit] to
            #end.'.format(time.ctime()))
            usr = input()  # break from schedule interval and manually call
            if usr.lower() =='exit':
                print('Goodbye.\n')
                thread.cancel()
                break
            else:
                thread.cancel()  # Stop the timer and cancel the execution of the
                #timer's action

    else:
        try:
            main(**options)
            # logging.shutdown()
        except Exception as e:
            logging.exception('Fatal error with Grouper function : {}'.format(e))
            raise
