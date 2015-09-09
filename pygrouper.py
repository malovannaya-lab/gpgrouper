# PyGrouper - Alex Saltzman
import re, csv, os, sys, itertools, time
import logging
import argparse
import threading
from collections import namedtuple, defaultdict
from math import ceil
from statistics import mean
# from sets import Set #python2 only
from configparser import SafeConfigParser, ParsingError, NoSectionError, NoOptionError
import pandas as pd
from pygrouper_subfuncts import *

try:
    from PIL import Image, ImageFont, ImageDraw
    imagetitle = True
except ImportError:
    imagetitle = False


def PeptidomeMatcher(usrdata, ref_dict):
    usrdata['metadatainfo'] = usrdata.apply(lambda x: genematcher(x['Sequence'], x['metadatainfo'], ref_dict), axis=1)

    return usrdata


def Grouper(usrfile, usrdata, exp_setup, FilterValues):
    global program_title
    global release_date

    print('Starting Grouper for exp file {}'.format(usrfile))
    logfilestr = '_'.join(str(x) for x in [exp_setup['EXPRecNo'],exp_setup['EXPRunNo'],exp_setup['EXPSearchNo'],\
                                           exp_setup['EXPTechRepNo'], str(exp_setup['EXPLabelType']) + '.log'])
    logfile = open(logfilestr, 'w+')  # append to previously opened log file
    logfile.write('{} | Starting {} for file : {}\n'.format(time.ctime(),program_title, usrfile))

    # -------------------- Populate gene, protein, homologene lists and counts -------------------- #
    usrdata['_data_tGeneList'], usrdata['_data_tProteinList'], usrdata['_data_GeneCount'], \
        usrdata['_data_ProteinCount'], usrdata['_data_tHIDList'], usrdata['_data_HIDCount'], \
        usrdata['_data_ProteinCapacity'] = list(zip(*usrdata.apply(lambda x: meta_extractor(x['metadatainfo']),axis=1)))
    # --------------------------------------------------------------------------------------------- #

    usrdata['_data_HID'], usrdata['_data_ProteinGI'], usrdata['_data_tTaxonIDList'], usrdata['_data_TaxonCount'] =\
                            '', '', '', ''  # potentially filled in later, fields will exist in database at least

    logfile.write('{} | Finished matching PSMs to {} refseq\n'.format(time.ctime(), usrdata.loc[0]['_data_TaxonID']))
    logging.info('Starting Grouper for exp number {}'.format(exp_setup['EXPRecNo']))
    nomatches = usrdata[usrdata['_data_GeneCount'] == 0].Sequence  # select all PSMs that didn't get a match to the refseq peptidome
    logfile.write('{} | Total identified PSMs : {}\n'.format(time.ctime(), \
                                                             len(usrdata[usrdata['_data_GeneCount'] != 0])))
    logfile.write('{} | Total unidentified PSMs : {}\n'.format(time.ctime(), len(nomatches)))

    for missing_seq in nomatches:
        # print('No match for sequence {} in {}').format(missing_seq,ursfile[0])  # No need to print this to screen
        logging.warning('No match for sequence {} in {}'.format(missing_seq, usrfile))
        # Store all of these sequences in the big log file, not per experiment.

    usrdata.drop('metadatainfo', axis=1, inplace=True)  # Don't need this anymore. But if we do, don't do this.
    resultout_dir = os.getcwd()  

    logfile.write('{} | Starting grouper.\n'.format(time.ctime()))
    usrdata['_data_PSM_IDG'] = usrdata.apply(lambda x: IDG_picker(x['IonScore'], x['q-Value']), axis=1) 
    #usrdata['Sequence'], usrdata['_data_SequenceModi'], usrdata['_data_SequenceModiCount'], usrdata['_data_LabelFLAG'] = list(zip(*usrdata.apply(lambda x : \
    #                                                                                                seq_modi(x['Sequence'], x['Modifications']), axis=1)))
    # We now move this to earlier, before we search against refseq (this will deal with the Xs
    # ==================== Gather all genes that all peptides map to ==================== #
    glstsplitter = usrdata['_data_tGeneList'].str.split(',').apply(pd.Series, 1).stack()
    glstsplitter.index = glstsplitter.index.droplevel(-1)  # get rid of multi-index
    glstsplitter.name = '_data_GeneID'  # give the series a name
    usrdata = usrdata.join(glstsplitter)  # usrdata gains column '_data_GeneID' from glstsplitter Series
    # =================================================================================== #

    logfile.write('{} | Starting peptide ranking.\n'.format(time.ctime()))
    usrdata = usrdata.sort(['Spectrum File', '_data_GeneID', 'Precursor Area', 'Sequence', 'Modifications',\
                            'Charge','_data_PSM_IDG','IonScore', 'PEP', 'q-Value'], 
                           ascending=[0, 1, 0, 1, 1, 1, 1, 0, 1, 1]) 
    usrdata.reset_index(inplace=True)
    usrdata.Modifications.fillna('', inplace=True)  # must do this to compare nans
    usrdata['Precursor Area'].fillna(0, inplace=True)  # must do this to compare nans
    grouped = usrdata.drop_duplicates(subset=['Spectrum File', '_data_GeneID', 'Sequence', 'Modifications', 'Charge']).groupby(['Spectrum File', '_data_GeneID'])  # each group belongs to 1 gene, now we can rank on a per-gene basis
    ranks = grouped.cumcount() + 1  # add 1 to start the peptide rank at 1, not 0
    ranks.name = '_data_PeptRank'
    usrdata = usrdata.join(ranks)
    usrdata['_data_PeptRank'] = usrdata['_data_PeptRank'].fillna(0)  # anyone who didn't get a rank gets a rank of 0
    logfile.write('{} | Peptide ranking complete.\n'.format(time.ctime()))
    print('{}: Peptide ranking complete for {}.'.format(datetime.now(), usrfile))
    logging.info('{}: Peptide ranking complete for {}.'.format(datetime.now(), usrfile))
    # =================================================================================== #

    usrdata['_data_AUC_nUseFLAG'], usrdata['_data_PSM_nUseFLAG'] =list(zip(*usrdata.apply(AUC_PSM_flagger, args=(FilterValues,), axis=1)))  # Flag good quality peptides
    pd.options.mode.chained_assignment = None  # default='warn'
    # Make the name for the peptide data table :
    usrdata_out = '_'.join(str(x) for x in [exp_setup['EXPRecNo'], exp_setup['EXPRunNo'],\
                                            exp_setup['EXPSearchNo'],  exp_setup['EXPTechRepNo'],\
                                            exp_setup['EXPLabelType'], 'data.tab'])
    # none/SILAC loop
    labeltypes = ['nolabel', 'heavy']  # right now only using nolabel, but in future this may grow
    for label in range(1):  # increase the range to go through more label types
        logfile.write('{} | Starting gene assignment for label {}.\n'.format(time.ctime(), labeltypes[label]))
        # =======================Select only peptides flagged  with good quality======================== #
        temp_df = usrdata[(usrdata['_data_LabelFLAG'] == label) & (usrdata['_data_AUC_nUseFLAG'] == 1) & (
            usrdata['_data_GeneCount'] > 0)]  # should keep WL's peptides too (if they have a mass)
        gene_cols = ['_e2g_EXPRecNo', '_e2g_EXPRunNo', '_e2g_EXPSearchNo', '_e2g_EXPLabelFLAG', '_e2g_AddedBy',
                     '_e2g_CreationTS', '_e2g_ModificationTS', '_e2g_GeneID', '_e2g_IDSet',
                     '_e2g_IDGroup','_e2g_IDGroup_u2g', '_e2g_GPGroup', '_e2g_GPGroups_All', '_e2g_PSMs',
                     '_e2g_PSMs_u2g', '_e2g_PeptidePrint', '_e2g_PeptideCount', '_e2g_PeptideCount_u2g',
                     '_e2g_PeptideCount_S', '_e2g_PeptideCount_S_u2g','_e2g_nGPArea_Sum_cgpAdj',
                     '_e2g_nGPArea_Sum_u2g', '_e2g_nGPArea_Sum_u2g_all', '_e2g_nGPArea_Sum_max',
                     '_e2g_nGPArea_Sum_dstrAdj', '_e2g_GeneCapacity', '_e2g_n_iBAQ_dstrAdj']  # cols of interest
        if args.genesets: gene_cols.append('_e2g_PeptideSet')  # this is sorta useless
        # =================================================================================== #
        #print(len(temp_df))
        if len(temp_df) > 0:  # only do if we actually have peptides selected
            genedata_out = '_'.join(str(x) for x in 
                                    [exp_setup['EXPRecNo'],exp_setup['EXPRunNo'],exp_setup['EXPSearchNo'],\
                                     exp_setup['EXPTechRepNo'], exp_setup['EXPLabelType'],\
                                     labeltypes[label], 'e2g.tab'])
            print('{}: Populating gene table for {}.'.format(datetime.now(), usrfile))
            # genes_df['_e2g_GeneID'] = Set(temp_df['_data_GeneID']) #only in 2.7
            genes_df = pd.DataFrame({'_e2g_GeneID': list(set(temp_df['_data_GeneID']))})
            genes_df['_e2g_GeneCapacity'] = genes_df.apply(lambda x: capacity_grabber(x['_e2g_GeneID'],\
                                                                                      temp_df), axis=1)
            genes_df['_e2g_EXPLabelFLAG'] = label
            for i, j in exp_setup.items(): genes_df[i] = j  # creates genes_df columns and populates them with 
                                                            # exp_setup values. Quick, easy, labeling

            genes_df['_e2g_PeptideSet'], genes_df['_e2g_PeptidePrint'], genes_df['_e2g_PeptideCount'], \
            genes_df['_e2g_PeptideCount_u2g'], genes_df['_e2g_PeptideCount_S'], genes_df['_e2g_PeptideCount_S_u2g'] = \
            list(zip(*genes_df.apply(PeptPrint1, args=(temp_df,), axis=1)))

            genes_df['_e2g_PSMs'], genes_df['_e2g_PSMs_u2g'], genes_df['_e2g_PSMs_S'],genes_df['_e2g_PSMs_S_u2g'] = \
            list(zip(*genes_df.apply(lambda x: e2g_PSM_helper(x['_e2g_GeneID'], \
                                                              temp_df, exp_setup['EXPTechRepNo']),axis=1)))

            print('{}: Calculating peak areas for {}.'.format(datetime.now(), usrfile))
            logging.info('{}: Calculating peak areas for {}.'.format(datetime.now(), usrfile))
            logfile.write('{} | Calculating peak areas.\n'.format(time.ctime()))

            genes_df['_e2g_nGPArea_Sum_max'], genes_df['_e2g_nGPArea_Sum_cgpAdj'], \
            genes_df['_e2g_nGPArea_Sum_u2g'], genes_df['_e2g_nGPArea_Sum_u2g_all']  = \
            list(zip(*genes_df.apply(AreaCalculator, args=(temp_df, exp_setup['EXPTechRepNo'], 
                                                           exp_setup['EXPQuantSource']), axis=1)))
                                                   # pandas gives a warning from this
                                                   # though I'm sure its fine

            print('{}: Calculating distributed area ratio for {}.'.format(datetime.now(), usrfile))
            logging.info('{}: Calculating distributed area ratio for {}.'.format(datetime.now(), usrfile))
            logfile.write('{} | Calculating distributed area ratio.\n'.format(time.ctime()))
            temp_df['_data_PrecursorArea_dstrAdj'] = temp_df.apply(AUC_distributor,args=(genes_df,
                                                                   exp_setup['EXPQuantSource'],), axis=1)
            if exp_setup['EXPQuantSource'] == 'AUC':
                normalize = 10 ** 9
            elif exp_setup['EXPQuantSource'] == 'Intensity':
                normalize = 10 ** 5

            print('{}: Assigning gene sets and groups for {}.'.format(datetime.now(), usrfile))
            logging.info('{}: Assigning gene sets and groups for {}.'.format(datetime.now(), usrfile))
            logfile.write('{} | Assigning gene sets and groups.\n'.format(time.ctime()))

            genes_df['_e2g_IDSet'], genes_df['_e2g_IDGroup'], genes_df['_e2g_IDGroup_u2g'] = list(
                zip(*genes_df.apply(GeneSetter, args=(genes_df, temp_df,), axis=1)))

            genes_df['_e2g_nGPArea_Sum_dstrAdj'] = genes_df.apply(gene_AUC_sum, args=(temp_df, normalize,), 
                                                                  axis=1)
            genes_df['_e2g_GeneCapacity'] = genes_df._e2g_GeneCapacity.astype('float')
            # print(genes_df._e2g_nGPArea_Sum_dstrAdj.dtype, genes_df._e2g_GeneCapacity.dtype)
            genes_df['_e2g_n_iBAQ_dstrAdj'] = genes_df._e2g_nGPArea_Sum_dstrAdj / genes_df._e2g_GeneCapacity
            genes_df['_e2g_GPGroup'] = ''
            genes_df.sort(columns=['_e2g_PSMs'], ascending=False, inplace=True)
            genes_df.index = list(range(0, len(genes_df)))
            last = 0
            for i in range(len(genes_df)):  # this is confusing but it works. The logic behind it is obvious, but its implementation is weird.
            # print last
                if genes_df.loc[i]['_e2g_IDSet'] != 3:
                    genes_df.loc[i]['_e2g_GPGroup'], lessthan = GPG_helper(genes_df.loc[i]['_e2g_IDSet'], \
                                                                          genes_df.loc[i]['_e2g_PeptideSet'], \
                                                                          genes_df, last)
                if type(genes_df.loc[i]['_e2g_GPGroup']) is int and not lessthan:
                    last = genes_df.loc[i]['_e2g_GPGroup']

            genes_df['_e2g_GPGroups_All'] = genes_df.apply(GPGall_helper, args=(genes_df,), axis=1)
            genes_df['_e2g_GPGroup'].replace(to_replace='', value=float('NaN'),
                                         inplace=True)  # can't sort int and strings, convert all strings to NaN
            genes_df.sort(columns=['_e2g_GPGroup'], ascending=True, inplace=True)
            genes_df.index = list(range(0, len(genes_df)))  # reset the index
            genes_df['_e2g_GPGroup'].fillna('', inplace=True)  # convert all NaN back to empty string
            # =============================================================#
            genes_df['_e2g_EXPRecNo'] = exp_setup['EXPRecNo']
            genes_df['_e2g_EXPRunNo'] = exp_setup['EXPRunNo']
            genes_df['_e2g_EXPSearchNo'] = exp_setup['EXPSearchNo']
            genes_df['_e2g_EXPTechRepNo'] = exp_setup['EXPTechRepNo']
            genes_df['_e2g_AddedBy'] = usrdata.loc[1]['_data_AddedBy']
            genes_df['_e2g_CreationTS'] = datetime.now().ctime()
            genes_df['_e2g_ModificationTS'] = datetime.now().ctime()
            genes_df.to_csv(genedata_out, columns=gene_cols, index=False, encoding='utf-8', sep='\t')
            logfile.write('{} | Export of genetable for labeltype {} completed.\n'.format(time.ctime(), 
                                                                                          labeltypes[label]))
    # -----------End of none/silac loop---------------------------------------- #  
    usrdata = pd.merge(usrdata, temp_df, how='outer')
    usrdata['_data_EXPRecNo'], usrdata['_data_EXPRunNo'], usrdata['_data_EXPSearchNo'], \
    usrdata['_data_EXPTechRepNo'] = exp_setup['EXPRecNo'], exp_setup['EXPRunNo'], exp_setup['EXPSearchNo'], \
    exp_setup['EXPTechRepNo']
    usrdata['_data_CreationTS'] = datetime.now().ctime()
    usrdata['_data_ModificationTS'] = datetime.now().ctime()
    usrdata['_data_tHIDList'] = ''  # will be populated later
    usrdata['_data_HIDCount'] = ''  # will be populated later
    data_cols = ['_data_EXPRecNo', '_data_EXPRunNo', '_data_EXPSearchNo', '_data_EXPTechRepNo', 'Sequence',
                 'PSM Ambiguity', 'Modifications', 'Activation Type', 'DeltaScore', 
                 'DeltaCn', 'Rank', 'Search Engine Rank', 'Precursor Area', 'QuanResultID', 'q-Value', 'PEP',
                 'Decoy Peptides Matched', 'Exp Value', 'Homology Threshold', 
                 'Identity High', 'Identity Middle', 'IonScore', 'Peptides Matched', '# Missed Cleavages',
                 'Isolation Interference [%]', 'Ion Inject Time [ms]', 'Intensity', 
                 'Charge', 'm/z [Da]', 'MH+ [Da]', 'Delta Mass [Da]', 'Delta Mass [PPM]', 'RT [min]',
                 'First Scan','Last Scan', 'MS Order', 'Matched Ions', 
                 'Total Ions', 'Spectrum File', 'Annotation', '_data_AddedBy', '_data_CreationTS',
                 '_data_ModificationTS', '_data_GeneID', '_data_tGeneList', '_data_ProteinGI', 
                 '_data_tProteinList','_data_ProteinCount', '_data_GeneCount', '_data_HID', 
                 '_data_tHIDList', '_data_HIDCount', '_data_TaxonID', '_data_tTaxonIDList', '_data_TaxonCount',
                 '_data_PSM_IDG', '_data_SequenceModi', 
                 '_data_SequenceModiCount', '_data_LabelFLAG', '_data_GeneID', '_data_PeptRank', 
                 '_data_AUC_nUseFLAG','_data_PSM_nUseFLAG', '_data_PrecursorArea_dstrAdj']
    #usrdata.to_csv(usrdata_out, columns=usrdata.columns, encoding='utf-8', sep='\t')
    #print(usrdata.columns.values)
    if not all(x in usrdata.columns.values for x in data_cols):
        print('Potential error, not all columns filled')
        data_cols = [x for x in data_cols if x in usrdata.columns.values]  # will still export successfully
    usrdata.to_csv(usrdata_out, columns=data_cols, index=False, encoding='utf-8', sep='\t')
    logfile.write('{} | Export of datatable completed.\n'.format(time.ctime()))
    logfile.write('Successful grouping of file completed.')
    logfile.close()
    print('Successful grouping of {} completed.\n' \
          .format('_'.join([str(exp_setup['EXPRecNo'])+'_'+str(exp_setup['EXPRunNo'])])))

def main(usrfiles=[], exp_setups=[], automated=False, usepeptidome=True):
    # ==========================Configuration Setup / Loading===================================#
    global program_title
    global release_date
    refs = {}
    if not os.path.isfile('py_config.ini'):
        input("No config file detected. Don't worry, we'll make one now\nPress [Enter] to continue")
        pysetup()
    elif args.setup:
        pysetup()
    parser = SafeConfigParser()
    parser.read('py_config.ini')
    pept_breakups = {}
    breakup_size = 4
    for taxon, location in parser.items('refseq locations'):
        refs[taxon] = {'loc': location, 'size': parser.getfloat('refseq file sizes',
                                                                taxon)}  # access and store preconfigured reference info
        if args.fullpeptread:
            pept_breakups[taxon] = [(0, int(refs[taxon]['size']))]
        else:
            pept_breakups[taxon] = [(0, int(ceil(refs[taxon]['size'] / breakup_size)))] * breakup_size
            # print pept_breakups
    # ==========================Configuration Setup / Loading===================================#

    if imagetitle:
        fancyprint(program_title, 12)  # ascii art
        #  fancyprint('Malovannaya lab',10)
    elif not imagetitle:
        print(program_title)  # not as exciting as ascii art

    print('\nrelease date: {}'.format(release_date))
    print('Python version ' + sys.version)
    print('Pandas version: ' + pd.__version__)

    if not automated: usr_name = input('Enter your name : ')

    FilterValues = {'Filter_IS': 7, 'Filter_qV': 0.05, 'Filter_PEP': 'all', 'Filter_IDG': 'all', 
                    'Filter_Z_min': 2,'Filter_Z_max': 4, 'Filter_Modi': 3}  # defaults
    Filter_Stamp = 'is{}_qv{}_pep{}_idg{}_z{}to{}mo{}'.format(FilterValues['Filter_IS'],
                                                              FilterValues['Filter_qV'] * 100,
                                                              FilterValues['Filter_PEP'],
                                                              FilterValues['Filter_IDG'], 
                                                              FilterValues['Filter_Z_min'],
                                                              FilterValues['Filter_Z_max'], 
                                                              FilterValues['Filter_Modi'])

    print('\nFilter values set to : {}'.format(Filter_Stamp))
    if not automated:
        try:
            explog = open('PyGrouper_grouped_exps.log', 'r+U')
        except IOError:
            explog = open('PyGrouper_grouped_exps.log', 'a+')
        grouped_exps = [value for value in re.findall(r'(\d+\S*\.txt)', ' '.join(explog.readlines()))]
        try:
            input('Press enter to continue, or Ctrl+C to modify the filter values.\n')
        except KeyboardInterrupt:
            print('\nFilter_PEP and Filter_IDG can be set to "all", but all other values must be numbers only.')
            dict_modifier(FilterValues, {'Filter_PEP': ['all'], 'Filter_IDG': ['all'], 'Filter_qV': []})

        logging.info('Filter stamp : {}'.format(Filter_Stamp))

        while True:
            try:
                exp_setup = {'taxonID': 9606, 'EXPTechRepNo': 1, 'EXPQuantSource': 'AUC', 'EXPRunNo': 1,\
                             'EXPSearchNo': 1, 'EXPLabelType': 'none'}
                usrfile_input = input('Enter a file to group or press Ctrl+C if done : ')
                if usrfile_input == 'forcequit':
                    logging.info('forcequit without selecting any files')
                    logging.shutdown()
                    sys.exit(1)
                if os.path.isfile(usrfile_input):  # check to see if file exists
                    proceed = True
                    if usrfile_input.strip() in grouped_exps:  # strip any whitespace to match correctly
                        while True:
                            proceed = input(
                                'Records show that {} has been grouped before. Are you sure you want to regroup (Y/n)? '\
                                .format(usrfile_input))
                            if 'Y' in proceed:
                                break
                            elif 'n' in proceed.lower():
                                proceed = False
                                break
                    if proceed:
                        try:
                            exp_setup['EXPRecNo'] = re.search('(\d{3,})',
                                                              usrfile_input).group()  # find first number of at least 3 digits
                        except AttributeError:
                            exprecno = input("Couldn't locate experimental record automatically, please enter it now. ")
                            exp_setup['EXPRecNo'] = exprecno
                        print('Experimental setup is : {}'.format(exp_setup))
                        try:
                            input('Press enter to accept values and continue, or Ctrl+C to modify the experimental values.\n')
                            exp_setups.append(exp_setup)
                        except KeyboardInterrupt:
                            print(
                                '\nEXPQuantSource can be set to AUC or Intensity, but all other values must be numbers only.')
                            dict_modifier(exp_setup, {'EXPQuantSource': ['AUC', 'Intensity']}, \
                                          exception_float=False)
                            exp_setups.append(exp_setup)

                        usrfiles.append(usrfile_input)
                else:
                    print('File {} not found.'.format(usrfile_input))
                    # morefile = raw_input('Do you have more files to group? ')
            except KeyboardInterrupt:
                if len(usrfiles) > 0:
                    print()
                    break
                else:
                    print('No files selected!')

    startTime = datetime.now()
    if FilterValues['Filter_PEP'] == 'all' or FilterValues['Filter_IDG'] == 'all':
        FilterValues['Filter_PEP'] = float('inf')

    print('\nStart at {}'.format(startTime))
    logging.info('Start at {}'.format(startTime))
    if args.fullpeptread: print('Running with fullpeptread option')
    usrdatas = []
    for usrfile in usrfiles:  # load all data
        usrdatas.append(pd.read_csv(usrfile, sep='\t'))
    for usrdata, exp_setup in zip(usrdatas, exp_setups):
        usrdata.rename(columns={'Annotated Sequence': 'Sequence', 'Area': 'Precursor Area',\
                                'Percolator q-Value': 'q-Value','Percolator PEP': 'PEP',\
                                'Ions Score': 'IonScore', 'DeltaM [ppm]': 'Delta Mass [PPM]',
                                'Deltam/z [Da]' : 'Delta Mass [Da]',
                                }, inplace=True)  # PD 2.0 column name changes
        if not automated:
            usrdata['_data_AddedBy'] = usr_name

        elif automated:
            usrdata['_data_AddedBy'] = exp_setup['AddedBy']

        usrdata['_data_TaxonID'] = exp_setup['taxonID']
        usrdata['_data_tGeneList'], usrdata['_data_tProteinList'], usrdata['_data_GeneCount'],\
            usrdata['_data_ProteinCount'], usrdata['_data_HomologeneID'], usrdata['_data_ProteinCapacity'], usrdata[
            'metadatainfo'] = '', '', 0, 0, '', '', ''

    # ==================== Load refseq and start matching peptides ==================== #
    print('{}: Loading refseq database.'.format(datetime.now()))
    RefseqInfo = namedtuple('RefseqInfo', 'taxonid, geneid, homologeneid, proteingi,genefraglen')
    for organism in refs.keys():
        if any(any(x['_data_TaxonID'].isin([int(organism)])) for x in
               usrdatas):  # check if we even need to read the peptidome for that organism
            # pept_reader = file_reader(refs[organism]['loc'])
            ref_reader = csv_reader(refs[organism]['loc'])
            print('Using peptidome {} '.format(refs[organism]))
            #print('Breakups : {}'.format(pept_breakups[organism]))
            #sys.exit(0)
            for breakup in pept_breakups[organism]:  # Read refseq in chunks, uses less memory
                prot = defaultdict(list)
                for k in range(*breakup):
                    try:
                        row = next(ref_reader)
                        fragments, fraglen = protease(row.fasta, minlen=7,
                                                      cutsites=['K', 'R'])  # should stop else clause here
                        # fragments = protease((linesplit[4].strip('\n'),minlen = 7, cutsites=['k','r'])
                        for fragment in fragments:
                            prot[fragment].append(RefseqInfo._make([row.taxonid, row.geneid, row.homologeneid, row.proteingi, \
                                                  fraglen]))

                    except StopIteration:  # breakups won't be exact since they are rounded up to ensure full file coverage ##should move this up and use else clause
                        break

                for usrdata, usrfile in zip(usrdatas, usrfiles):
                    #print(usrdata.loc[0]['_data_TaxonID'])
                    if usrdata.loc[0]['_data_TaxonID'] == int(organism):  # check to see if the inputdata 
                                                                         #taxonID matches with the proteome
                        usrdata['Sequence'], usrdata['_data_SequenceModi'], usrdata['_data_SequenceModiCount'], usrdata['_data_LabelFLAG'] = \
                        list(zip(*usrdata.apply(lambda x : seq_modi(x['Sequence'], x['Modifications']), axis=1)))
                        # print 'Matching for {}'.format(usrfile)
                        PeptidomeMatcher(usrdata, prot)  # call matcher
                del prot  # free up memory

    print('{}: Finished matching peptides to genes.'.format(datetime.now()))
    logging.info('{}: Finished matching peptides to genes.'.format(datetime.now()))
    failed_exps = []
    for usrfile, udata, esetup in zip(usrfiles, usrdatas, exp_setups):
        # print usrfile, esetup
        try:
            Grouper(usrfile, udata, esetup, FilterValues)
        except Exception as e:  # catch and store all exceptions, won't crash whole program at least
            failed_exps.append((usrfile, e))
            print('Failure for file of experiment {}.\nThe reason is : {}'.format(esetup['EXPRecNo'], e))
            logging.warn('Failure for file of experiment {}.\nThe reason is : {}'.format(esetup['EXPRecNo'], e))
            raise  # usually don't need to raise, will kill the script. Re-enable if need to debug and find where errors are
    print('Time taken : {}'.format(datetime.now() - startTime))
    logging.info('Time taken : {}.\n\n'.format(datetime.now() - startTime))
    if automated:
        return failed_exps
    else:
        if failed_exps:  # list of failed experiments, is empty if no failures
            for failed in failed_exps:  # failed is a tuple with the datafile and then the reason for failure
                usrfiles.remove(failed[0])  # remove all of the failed experiments so they won't go in log
                failedlog.write('{} : failed grouping {}, reason : {}\n'.format(datetime.now(), *failed))
        if usrfile:
            for usrfile in usrfiles:
                explog.write('{} : grouped experiment file {} \n'.format(datetime.now(), usrfile))


def expchecker(options):
    '''Checks for new experimental data based on 
    '''
    global program_title
    # exp_setup = {'taxonID' : 9606, 'EXPTechRepNo' : 1, 'EXPQuantSource' : 'AUC', 'EXPRunNo' : 1, 'EXPSearchNo': 9}
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
    # grouped_exps = [int(value) for value in re.findall(r'(\d+)',' '.join(explog.readlines()))]
    grouped_exps = [value for value in re.findall(r'(\d+\S*\.txt)', ' '.join(
        explog.readlines()))]  # starts with at least 1 digit, has unspecified amount of non-whitespace, and ends in .tab
    failed_exps = [value for value in re.findall(r'(\d+\S*\.txt)', ' '.join(failedlog.readlines()))]
    # failed_exps = [int(value) for value in re.findall(r'(\d+)',' '.join(failedlog.readlines()))]
    validfiles = [f for f in os.listdir(os.getcwd()) if '.txt' in f]
    ispec_export = '4PyGrouper_ExpRunDump.xlsx'
    experiment_file = pd.read_excel(ispec_export)
    # gen = experiment_file.iterrows()
    exp_setups, usrfiles = [], []
    required_fields = ['_exprun_EXPRecNo','_exprun_EXPRunNo']
    experiment_file.fillna(0, inplace=True)  # Fill NAs with 0 so we can convert everything to int
    for field in required_fields:
        experiment_file[field] = experiment_file[field].astype('int')  # get rid of decimals
    for row in experiment_file.iterrows():  # row is a tuple, [0] is index, [1] is actual row
        if not any(row[1][required_fields].isin([0])):  # ensure all required info is present
            record = str(row[1]._exprun_EXPRecNo) + '_' + str(row[1]._exprun_EXPRunNo)
            try:
                if all(record not in exp for exp in grouped_exps) and all(
                            record not in fexp for fexp in failed_exps):  # then we group
                    exp_setup = {'EXPRecNo': int(row[1]._exprun_EXPRecNo), 'EXPRunNo': 
                                 int(row[1]._exprun_EXPRunNo),
                                 'EXPSearchNo': int(row[1]._exprun_EXPSearchNo),
                                 'taxonID': int(row[1]._exprun_TaxonID),
                                 'EXPQuantSource': str(row[1].exprun_Search_QuantSource),
                                 'AddedBy': str(row[1]._exprun_AddedBy), 'EXPTechRepNo': 
                                 int(row[1].exprun_nTechRepeats),
                                 'EXPLabelType': str(row[1].exprun_LabelType)}
                    passed = True  # flag for if we should group an experiment
                else:
                    passed = False
            except ValueError:
                # print 'Experimental information not complete for experiment number 
                passed = False
                pass
        
            if passed:
                usrfilelst = [x for x in validfiles if record in x]  # match exp record and run numbers
                if len(usrfilelst) == 1:  # ensure we just have 1 file
                    usrfile = usrfilelst[0]
                elif len(usrfilelst) > 1:
                    usrfilelst = [x for x in usrfilelst if 'all' in x]  # find the file with 'all' in it,
                                                             # this is based on how we name files for MudPIT
                    if len(usrfilelst) == 1:
                        usrfile = usrfilelst[0]
                    else:  # can't find the correct file if we get to here
                        print('Warning, more than one match for EXPRecNo {} detected, skipping...'.format(
                            exp_setup['EXPRecNo']))
                        passed = False
                elif len(usrfilelst) == 0:  # no file found
                # print 'Warning, no file found for EXPRecNo {}.'.format(exp_setup['EXPRecNo'])
                    passed = False
            if passed and len(exp_setups) <= 25:  # cap on how may files to group at once, untested
                exp_setups.append(exp_setup)
                usrfiles.append(usrfile)
                passed = False  # reset flag

    if len(usrfiles) > 0:  # then we group
        for usrfile, exp_setup in zip(usrfiles, exp_setups):
            print('Found experiment {} from datafile {}.'.format(exp_setup['EXPRecNo'], usrfile))
            logfilestr = '_'.join(
                str(x) for x in [exp_setup['EXPRecNo'], exp_setup['EXPRunNo'], exp_setup['EXPSearchNo'],
                                 exp_setup['EXPTechRepNo'], exp_setup['EXPLabelType'] + '.log'])
            logfile = open(logfilestr, 'w')
            logfile.write('{} | starting {} for {} with taxonID: {}\n'.format(time.ctime(), program_title, 
                                                                              exp_setup['EXPRecNo'],      
                                                                              exp_setup['taxonID']))
            logfile.close()

        try:
            logging.info('{}: Initiating {}'.format(datetime.now(), program_title))
            failed_exps = main(usrfiles, exp_setups, **options)
            if failed_exps:  # list of failed experiments, is empty if no failures
                for failed in failed_exps:  # failed is a tuple with the datafile and then the reason for failure
                    usrfiles.remove(failed[0])  # remove all of the failed experiments so they won't go in log
                    failedlog.write('{} : failed grouping {}, reason : {}\n'.format(datetime.now(), *failed))
            if usrfile:
                for usrfile in usrfiles:
                    explog.write('{} : grouped experiment file {} \n'.format(datetime.now(), usrfile))
        except Exception as e:  # catch all exceptions
            logging.exception('Fatal error with Grouper function : {}'.format(e))
            raise  # crash and burn
    failedlog.close()
    explog.close()
    #print('{} : Sleeping... Press Enter to wake or Ctrl+C to end.'.format(time.ctime()))


def schedule(INTERVAL, options):
    expchecker(options)
    global thread
    thread = threading.Timer(INTERVAL, schedule, [INTERVAL, options])
    print('{} : Sleeping... Press Enter to wake or [exit] to end.'.format(time.ctime()))
    thread.start()


if __name__ == '__main__':
    program_title = 'PyGrouper v0.1.004'
    release_date = '11 August 2015'
    parser = argparse.ArgumentParser()
    parser.add_argument('-fpr', '--fullpeptread',
                        action='store_true',
                        help='Default FALSE. Load the peptidome all at once, minor speed improvement at the expense of memory. Use at your own risk.')
    parser.add_argument('-gs', '--genesets',
                        action='store_true',
                        help='Optional inclusion of peptide sets for each gene. Legacy and not useful.')
    parser.add_argument('-s', '--setup',
                        action='store_true',
                        help='Run setup wizard for PyGrouper.\nNot necessary unless adding a new refseq or need to change refseq location.')
    parser.add_argument('-auto', '--automated',
                        action='store_true', help='Automated run of PyGrouper. Note, requires experiment dump from iSPEC to be set up seperately.\nIf you are not sure if this is correctly set up, automation will probably not work.')
    args = parser.parse_args()
    logfilename = program_title.replace(' ', '_') + '.log'
    logging.basicConfig(filename=logfilename, level=logging.DEBUG)
    logging.info('{}: Initiating {}'.format(datetime.now(), program_title))
    options = {}
    if args.automated:
        options['automated'] = True
        while True:
            
            INTERVAL = 60 * 60
            schedule(INTERVAL, options)
            #print('{} : Sleeping... Press Enter to wake or [exit] to end.'.format(time.ctime()))
            usr = input()  # break from schedule interval and manually call
            if usr.lower() =='exit':
                print('Goodbye.\n')
                thread.cancel()
                break
            else:
                thread.cancel()  # Stop the timer and cancel the execution of the timer's action

    else:
        try:
            main(**options)
            # logging.shutdown()
        except Exception as e:
            logging.exception('Fatal error with Grouper function : {}'.format(e))
            raise
