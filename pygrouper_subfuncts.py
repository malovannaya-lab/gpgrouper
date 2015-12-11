#===============================================================================#
import re, sys
import os
import csv
import numpy as np
from datetime import datetime
from collections import namedtuple, defaultdict, OrderedDict
from configparser import ConfigParser
from statistics import mean
import logging
import json
from collections import defaultdict
try : 
    from PIL import Image, ImageFont, ImageDraw
    #imagetitle = True
except ImportError: pass
    #imagetitle = False

RefseqInfo = namedtuple('RefseqInfo',
                        'taxonid, geneid, homologeneid,proteingi,genefraglen')

def column_constructor(columndict=None):
    print(columndict)
    if not columndict:
        columndict = defaultdict(list)
    columns = OrderedDict([
        ('Sequence','the peptide sequence'),
        ('Modifications', 'post translational modifications'),
        ('PSM_Ambiguity', ''),
        ('ActivationType', 'Activation method for MS2'),
        ('DeltaScore', ''),
        ('DeltaCn', ''),
        ('Rank', ''),
        ('SearchEngineRank', ''),
        ('PrecursorArea','MS1 integrated peptide area'),
        ('Intensity','peptide intensity'),
        ('q_value','FDR adjusted q-value'),
        ('PEP','Posterior Error Probability'),
        ('IonScore', 'Mascot Ion Score'),
        ('MissedCleavages', ''),
        ('IsolationInterference', ''),
        ('IonInjectTime', ''),
        ('Charge', 'Ion charge'),
        ('m/z [Da]', ''),
        ('MH+ [Da]', ''),
        ('MatchedIons', ''),
        ('TotalIons', ''),
        ('Delta Mass [Da]',''),
        ('Delta Mass [PPM]',''),
        ('RT [min]', ''),
        ('MS Order', ''),
        ('SpectrumFile', ''),
        
        ('psm_EXPRecNo', 'Experiment record number'),
        ('psm_EXPRunNo', 'Experiment run number'),
        ('psm_EXPSearchNo', 'Experiment search number'),
        ('psm_EXPTechRepNo', 'Experiment tech rep number'),
        ('psm_AddedBy', 'Person who added the experiment'),
        ('psm_CreationTS', 'experiment creation timestamp'),
        ('psm_ModificationTS', 'experiment modification timestamp'),
        ('psm_GeneID', 'geneID for a given peptide'),
        ('psm_GeneList', 'all genes that a given peptide maps to'),
        ('psm_GeneCount', 'count of the number of genes for a given peptide'),
        ('psm_ProteinGI', 'ProteinGI for a given peptide'),
        ('psm_ProteinList', 'all proteinGIs for a given peptide'),
        ('psm_ProteinCount', 'number of proteins for a given peptide'),
        ('psm_HID', 'homologene ID for a given peptide geneID'),
        ('psm_HID_list','list of all homologenes for a given peptide'),
        ('psm_HID_count', 'count of all homologenes for a given peptide'),
        ('psm_TaxonID','taxonID for a given peptide'),
        ('psm_TaxonIDlist', 'list of all taxonIDs for a given peptide,'\
          ' useful when searching against multiple species'),
        ('psm_TaxonCount', 'count of the number of taxons for a given peptide'),
        ('psm_PSM_IDG', 'PSM IDG for a given peptide'),
        ('psm_SequenceModi', 'string of sequence with modifications'),
        ('psm_SequenceModiCount', 'count of total modifications'),
        ('psm_ModiCount', 'number of modifications for a given peptide'),
        ('psm_LabelFLAG', 'label number for a given peptide'),
        ('psm_PeptideRank', 'rank for a given peptide within its geneID'),
        ('psm_AUC_useflag', 'boolean value for use of peptide for calculating area'),
        ('psm_PSM_useflag', 'boolean value for use of peptide for psm'),
        ('psm_Area_taxonAdj','taxon adjusted area'),
        ('psm_PrecursorArea_dstrAdj','redistributed peptide precursor area'),
        
        ('gene_EXPRecNo', 'Experiment Record Number'),
        ('gene_EXPRunNo', 'Experiment Run Number'),
        ('gene_EXPSearchNo', 'Experiment Search Number'),
        ('gene_EXPLabelFLAG', 'label type, eg. None, TNT, etc'),
        ('gene_AddedBy', 'Person who added the experiment'),
        ('gene_CreationTS', 'experiment creation timestamp'),
        ('gene_ModificationTS', 'experiment modification timestamp'),
        ('gene_GeneID', 'gene product GeneID'),
        ('gene_IDSet', 'IDSet'),
        ('gene_IDGroup','IDGroup'),
        ('gene_IDGroup_u2g','IDGroup for unique to gene peptide'),
        ('gene_GPGroup', 'GPGroup'),
        ('gene_GPGroups_All','GPGroups All'),
        ('gene_PSMs', 'number of psms for that gene product'),
        ('gene_PSMs_u2g', 'number of unique psms for that gene product'),
        ('gene_PeptidePrint', 'Print of the peptides identified for a given gene product'),
        ('gene_PeptideCount', 'count of peptides for a given gene product'),
        ('gene_PeptideCount_u2g', 'count of unique peptides for a given gene product'),
        ('gene_PeptideCount_S', 'count of strict peptides for a given gene product'),
        ('gene_PeptideCount_S_u2g','count of strict, unique peptides for a given gene product'),

        ('gene_GeneArea_gpcAdj','gene product area adjusted for gene product counts'),
        ('gene_GeneArea_Sum_cgpAdj', ''),
        ('gene_GeneArea_gpcAdj_u2g','gene product area adjusted for gene product counts, u2g'),
        ('gene_GeneArea_gpcAdj_u2g_all','description goes here'),
        ('gene_GeneArea_gpcAdj_max','maximum potential gene product area'),
        ('gene_GeneArea_dstrAdj', 'gene product area after peptide area distribution'),
        ('gene_GeneCapacity', 'number of tryptic peptides for a given gene product'),
        ('gene_iBAQ', 'dstrAdj gene area divided by the gene capacity'),

    ])

                
    existing_values = [item for sublist in
                       columndict.values() for
                       item in sublist]
    proceed = True
    skipall = False
    for col in columns:
        if not proceed:
            break
        print()
        if col not in existing_values:
            columndict[col].append(col)
            existing_values.append(col)
        while True and not skipall:
            print('Enter an alias for '\
                  'the column name {} or '\
                  'type [Done] if finished.\n'\
                  'Type [description] to get a description '\
                  'of a given column.'.
                  format(col))
            print('Type done without entering anything to use the '\
                  'default name.')
            print('')
            print('The current values are : {}'.format(
                ', '.join(columndict[col])))
            newValue = input('>>> ')
            if newValue.strip().lower() == 'exit':
                proceed = False
                break
            elif newValue.strip().lower() == 'done':
                break
            elif newValue.strip().lower() == 'description':
                print('{} : {}\n'.format(col, columns[col]))
                input('[OK')
            elif newValue == 'skipall':
                skipall = True
                break
            elif newValue and newValue not in existing_values:
                columndict[col].append(newValue)
                existing_values.append(newValue)
            elif newValue in existing_values:
                print('Alias already exists! '\
                        'Please choose another.\n')
                input('[OK]')
            elif not newValue:
                print('Invalid entry!\n')

    return columndict

def column_to_ini(BASE=''):
    parser = ConfigParser()
    parser.optionxform = str
    if os.path.isfile('py_config.ini'):
        parser.read('py_config.ini')
    elif not os.path.isfile('py_config.ini'):
        print('No config file found.')
        pysetup()
        parser.read('py_config.ini')
    if 'column names' not in parser.sections():
        parser.add_section('column names')
    column_aliases = defaultdict(list)  # load previous 
    for key in parser['column names'].keys():
        value = parser.get('column names',key)
        aliases = list(filter(None,(x.strip() for
                                    x in value.splitlines())))
        column_aliases[key] = aliases
    #column_aliases = dict(parser['column names'])
    column_aliases = column_constructor(columndict=column_aliases)
    print(column_aliases)
    for key in column_aliases:
        value = '\n'+'\n'.join(column_aliases[key])
        parser.set('column names', key, value)
    parser.write(open('py_config.ini','w'))    
    return parser
    #return column_aliases

def column_identifier(df, aliases):
    column_names = dict()
    for col in aliases:
        for alias in aliases[col]:
            name = [dfcolumn for dfcolumn in df.columns if dfcolumn==alias]
            if len(name)==1:
                column_names[col] = name[0]
                break
    return column_names
       
        

def json_column_constructor(BASE='',append=False ):
    '''Function for reading (and creating if necessary) a json
       file for mapping custom column names to a universal name set'''
    jsonFile = os.path.join(Base,'column_headers.json')
    to_json = False
    if os.path.isfile(jsonFile):
         columndict = json.load(jsonFile)
         if append:
            columndict = column_constructor(columndict)
            to_json = True
    elif not os.path.isfile(jsonFile):
       columndict = column_constructor()
       to_json = True
    if to_json:
        json.dump(os.path.join(BASE,columndict),
                   indent=4, sort_keys=True)
    columndict = json.load(jsonFile)
    return
    
def bufcount(filename):
    '''fast way to count lines in a file'''
    f = open(filename, mode = 'rt', encoding = 'latin1')                  
    lines = 0
    buf_size = 2048
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines+1  # plus 1 since I use csv_reader, which ignores the header

    
def pysetup():
    
    parser = ConfigParser()
    
    cancel = False
    sections = []
    if os.path.isfile('py_config.ini'):
        parser.read('py_config.ini')
        sections = parser.sections()
        #while True:
            #del_confirm = input('A config file is already detected? '\
            #   'Are you sure you want to delete and re-run the setup(Y/n)? ')
            #if del_confirm == 'Y' :
                #    os.remove('py_config.ini')
                #    break
            #elif 'n' in del_confirm.lower() : cancel = True
    if not cancel:    
        for section in ['refseq locations', 'refseq file sizes']:
            if section not in sections:
                parser.add_section(section)
        print('Welcome to PyGrouper setup wizard. First we need to get a '\
              'list of all taxons you wish to be able to search for.')
        print('Note this requires that you have refseq data for each taxon'\
              ' you list.')
        taxons = set(parser['refseq locations'].keys())
        while True:
            try:
                if len(taxons) > 0 :
                    print('\nCurrently added taxons : '\
                             ' {}\n'.format(' '.join(taxons)))
                newtaxon = input('Enter a taxon ID to add (human is 9606) '\
                                 'or Ctrl+C if done : ')
                taxons.add(newtaxon)
            except KeyboardInterrupt:
                break
        for taxon in taxons:
            print()
            cancel = False
            if parser['refseq locations'].get(taxon):
                print('Taxon {} already has a location assigned : {}'.format(
                    taxon, parser.get('refseq locations',taxon)))
                reassign = input('Would you like to reassign? [Y/n] ')
                if reassign != 'Y':
                    cancel = True
            while True and not cancel:
                taxon_loc = input('\nEnter location of the refseq file'\
                                  ' for taxon {} : '.format(taxon))
                if os.path.isfile(taxon_loc):
                    print('Setting {} refseq location to {}'.format(taxon,
                                                                    taxon_loc))
                    parser.set('refseq locations',taxon,taxon_loc)
                    print('Calculating...')
                    parser.set('refseq file sizes',taxon,
                               str(bufcount(taxon_loc)))
                    print('Done!\n')
                    break
                else : 
                    print('\nFile not found!')
                    if input("Type 'pass' to skip taxon {} from list "\
                             "or [Enter] to try again : ".format(
                                  taxon)).lower() == 'pass':
                        break
        parser.write(open('py_config.ini','w'))
        print('Configuration of PyGrouper is complete.\n')
    elif cancel:
         print('Canceling..\n')
   

def csv_reader(inputfile,sep='\t'):
    FastaRecord = namedtuple(
         'FastaRecord','taxonid, homologeneid, geneid, proteingi, fasta')
    with open(inputfile, mode ='rU') as f:
        reader = csv.DictReader(f,delimiter = sep)
        for line in reader:
            record = FastaRecord._make([line['faa_TaxonID'],
                                        line['faa_HomologeneID'],
                                        line['faa_GeneID'],
                                        line['faa_ProteinGI'],
                                        line['faa_FASTA']])
            yield record


def protease(seq,minlen = 0,cutsites = [], exceptions = []):
    frags = []
    chop = ''
    while seq:
        cuts = [seq.find(x) for x in cutsites]
        if len(cuts) > 0 :
            if min(cuts) >= 0 : cut = min(cuts) + 1
            else : cut = max(cuts) + 1
        else : cut = 0
        chop += seq[:cut]
        seq = seq[cut:] 
        #print(chop, cut, seq)
        if cut == 0 or len(seq) == 0: 
            #print(chop)
            if len(seq) != 0 :frags.append(seq)
            elif len(seq) == 0 and chop : frags.append(chop) #special case for
                                   #if last amino acid is one of the cut sites
            break #no more cuts, wrap up and go home
        if seq[0] not in exceptions:
            frags.append(chop)
            #print(chop, seq)  # for debugging
            chop = ''
        
    frags_1, frags_2  = [], [] 
    for i in range(len(frags)-1):
        frags_1.append(frags[i]+frags[i+1])
    for i in range(len(frags)-2):
        frags_2.append(frags[i]+frags[i+1]+frags[i+2])
    #for i in range(len(frags)-3):
        #frags_2.append(frags[i]+frags[i+1]+frags[i+2]+frags[i+3])    
    #print frags, frags_1, frags_2
    [a.append(a[0][1:]) for a in [frags, frags_1, frags_2,] if len(a) > 0]
    # chop off each N-term
    frags = [x for x in frags if len(x) >= minlen]
    fragments = [x for x in frags + frags_1 + frags_2 if len(x) >= minlen]  
    return fragments, len(frags) #len(frags) is no miscuts
    
def genematcher(seq, metadata, prot):
    seq = seq.upper()
    #print(type(metadata))  # for debugging
    if not isinstance(metadata,list):
        metadata = []
    if seq in prot:
        for data in prot[seq]:
            metadata.append(data)
    return metadata

def genelist_extractor(metadata):
     ''' at the peptide level'''
     taxonids = set()
     homologeneids = set()
     proteingis = set()
     genefraglens = []

     genelist = set()
     if metadata:
          for data in metadata:
             genelist.add(data.geneid)
             taxonids.add(data.taxonid)
             homologeneids.add(data.homologeneid)
             proteingis.add(data.proteingi)
             #gene_metadata[data.geneid].append((data.taxonid, data.homologeneid,
             #                             data.proteingi, data.genefraglen))

     #for gene in genes:
     #     if gene not in genelist:  # we split based on genelist, need unique set
     #          genelist.append(gene)

     return (','.join(genelist), len(genelist), ','.join(taxonids),
             len(taxonids), ','.join(proteingis),
             len(proteingis))

def meta_extractor(geneid, genedict):  # need to take in whole tempdf.
                                       # not at peptide level anymore
     ''' at the gene level 
     taxonid, homologeneid, proteingi, genefraglen'''

     homologeneids = set()
     proteingis = set()
     genefraglens = []
     for data in genedict[geneid]:

          homologeneids.add(data[1])
          proteingis.add(data[2])
          genefraglens.append(data[3])

     if not genefraglens:
          genefraglens.append(0)  # if nothing gets populated

     return (','.join(homologeneids),
             len(homologeneids), ','.join(proteingis), len(proteingis),
             mean(genefraglens))
             
     

def meta_extractor_old(metadata):
    taxonids = []
    genes = defaultdict(list)
    genelst = []
    homologenes = []
    proteingis = []
    genefraglen = []
    frags = ['0'] # will be overwritten if there are genes    
    if metadata:
        for data in metadata:
            genes[data.geneid].append((data.taxonid, data.homologeneid,
                                       data.proteingi, data.genefraglen))
        #print(genes)
        #sys.exit(1)
            #genelst.append(data.geneid)
            #taxonids.append(data.taxonid)
            #homologenes.append(data.homologeneid)
            #proteingis.append(data.proteingi)
            #genefraglen.append(str(data.genefraglen))

        for gene in genes:
            if gene not in genelst:
                genelst.append(gene)
                frags = []
                for data in genes[gene]:
                #if data[0] not in taxonids:
                    taxonids.append(data[0])  # next step
                    frags.append(str(data[3]))
                #if data[1] not in homologenes:

                    homologenes.append(data[1])
                #if data[2] not in proteingis:
                    proteingis.append(data[2])

                #if len(data[3]) > 0: 
                     
                #else:
                #     frags.append('0')

        #if len(frags) > 0 : genefraglen.append(str(mean(frags)))
        #else : genefraglen.append('0')
    #print(genefraglen)
    return ','.join(genelst), ','.join(proteingis), len(genelst),\
         len(proteingis), ','.join(homologenes), len(homologenes),\
         ','.join(frags), ','.join(taxonids), len(taxonids)             

def dict_modifier(d,exceptions = {},exception_float = True):
    """Dictionary modifier function for the user. Takes in dictionary (d), as 
    well as all dictionary keys that do not have to be floats (exceptions),
    the valid strings that the exception key values can hold (valid_exceptions), 
    as well as a boolean value for if the exception key values can be numbers
    """
    while True:
        print('\nThe current settings are : ')
        for key,value in list(d.items()): print('{} : {}'.format(key,value))
        try:
            
            filter_change = input('Type a value to change (case sensitive),'\
                                  ' or hit Ctrl+C if done. ')
            try :
                new_value = input('The current value is set to {}. Type a'\
                                  ' new value. : '.format(d[filter_change]))
                if filter_change not in list(exceptions.keys()):
                    try : 
                        new_value_float = int(new_value)
                        logging.info('{} changed from {} to {}'.format(
                             filter_change, d[filter_change], new_value_float))
                        d[filter_change] = new_value_float   
                    except ValueError:
                        print('\nInvalid entry. New value must be an integer.\n')
                elif filter_change in list(exceptions.keys()):
                    if new_value in exceptions[filter_change]:
                        logging.info('{} changed from {} to {}'.format(
                             filter_change, d[filter_change], new_value))
                        d[filter_change] = new_value
                        
                    elif exception_float:
                        try : 
                            new_value_float = float(new_value)
                            logging.info('{} changed from {} to {}'.format(
                                 filter_change, d[filter_change],
                                 new_value_float))
                            d[filter_change] = float(new_value_float)
                        except ValueError:
                            print('\nInvalid entry. New value must be a '\
                                  'number or in {}.\n'.format(valid_exceptions))
                    else : print('\nInvalid entry. New value must be in'\
                                 ' {}.\n'.format(valid_exceptions))
            except KeyError :
                print('\nInvalid filter value.\n')
                
        except KeyboardInterrupt:
            print('')
            break
    return d   

def nan_popper(lst):
    return [x for x in lst if not np.isnan(x)]

    
def fancyprint(ShowText,string_size=12): 
    try :
        font = ImageFont.truetype('arialbd.ttf', string_size) #load the font
        font_import  = True
    except IOError:
        font_import = False
        print(ShowText)
    if font_import:    
        size = font.getsize(ShowText)  #calc the size of text in pixels
        image = Image.new('1', size, 1)  #create a b/w image
        draw = ImageDraw.Draw(image)
        draw.text((0, 0), ShowText, font=font) #render the text to the bitmap
        for rownum in range(size[1]): 
        # scan the bitmap:
        # print ' ' for black pixel and 
        # print '#' for white one
            line = []
            for colnum in range(size[0]):
                if image.getpixel((colnum, rownum)): line.append(' '),
                else: line.append('#'),
            print(''.join(line))      

#def AUC_PSM_flagger(df,Filter_IS, Filter_qV,Filter_PEP, Filter_IDG, Filter_Z_min, Filter_Z_max, Filter_Modi):
def AUC_PSM_flagger(df,d):
    if d['Filter_PEP'] =='all' : d['Filter_PEP'] = float('inf')
    if d['Filter_IDG'] =='all' : d['Filter_IDG'] = float('inf')
    #if df['IonScore'] == '' and df['q-Value'] == '': 
    if df['Charge'] < d['Filter_Z_min'] or df['Charge'] > d['Filter_Z_max'] :
         # if charge is empty (nan), this will not be true
        AUC_flag = 0
        PSM_flag = 0

    elif df['psm_SequenceModiCount'] > d['Filter_Modi'] : 
        AUC_flag = 0
        PSM_flag = 0

    elif np.isnan(df['IonScore']) and np.isnan(df['q_value']): 
        AUC_flag = 1 #retains WLs Q2up assignments
        PSM_flag = 0

    elif df['IonScore'] < d['Filter_IS'] : 
        AUC_flag = 0
        PSM_flag = 0
    elif df['q_value'] > d['Filter_qV'] : 
        AUC_flag = 0
        PSM_flag = 0
    elif df['PEP'] > d['Filter_PEP'] : 
        AUC_flag = 0
        PSM_flag = 0
    elif df['psm_PSM_IDG'] > d['Filter_IDG'] : 
        AUC_flag = 0
        PSM_flag = 0
    elif df['psm_PeptideRank'] == 0 : 
        AUC_flag = 0 # omit multiPSM peak areas
        PSM_flag = 1  # Change from 0 in grouper
    else: 
        AUC_flag = 1
        PSM_flag = 1
    return AUC_flag, PSM_flag
    

def IDG_picker(IonScore, qvalue):

    if (IonScore>= 30 and qvalue <= .01): IDGout = 1
    elif (IonScore >= 30 and qvalue <= .05): IDGout = 2
    elif (IonScore >= 20 and qvalue <= .01): IDGout = 3
    elif (IonScore >= 20 and qvalue <= .05): IDGout = 4
    elif (IonScore >= 10 and qvalue <= .01): IDGout = 5       
    elif (IonScore >= 10 and qvalue <= .05): IDGout = 6
    elif (IonScore >= 0 and qvalue <= .01): IDGout = 7
    elif (IonScore >= 0 and qvalue <= .05): IDGout = 8
    else: IDGout  = 9
    return IDGout
 
        
def seq_modi(sequence, modifications):

    '''
    function to output a modified string of sequence that includes all 
    modifications at the appropriate place.
    Modified amino acids are indicated in the input 'sequence' by a lowercase 
    letter.
    A dictionary of potential modification keys and their corresponding 
    substitutions in the primary sequence is provided;
    all other modifications should be ignored.
    Sometimes N-terminal modifications are present, this should be considered 
    the same as a modification at the first amino acid.
    '''
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    amino_acids = ['('+x+')' for x in amino_acids]
    modtext = {'(DeStreak)' : '(des)', '(Deamidated)' : '(dam)',
               '(Carbamidomethyl)' : '(car)', '(Oxidation)' : '(oxi)',
               '(Phospho)' : '(pho)', '(Prot)(Acetyl)':'(ace)',
               '(Acetyl)': '(ace)', '(GlyGly)' : '(gg)', '(Label:13C(6))' :
               '(lab)', '(Label:13C(6)+GlyGly)' : '(labgg)'}
    seqmodi = ''
    seqlist = list(sequence)
    modi_len = 0  # default length is zero, can change
    label = 0  # default label = 0, can change
    if any(c for c in sequence if (c.islower() or c == 'X')):
         # check if any modifications to deal with
        if 'Label' in modifications:
             #routine for grabbing everything enclosed in parenthesis, including
             #nested. Easier to write than regex when dealing with nested
            label = 1  # for indication of label on peptide
            openb, modi, modkeys = 0, '', []
            for e in modifications:
                if e == '(':
                    openb += 1
                elif e == ')':
                    openb  += -1
                if openb:
                    modi +=e
                elif not openb:
                    if modi:
                        modi += e  #add the last parenthesis
                        modkeys.append(modi)
                    modi = ''
        else:
            modkeys = re.findall(r'(\([^\)]+\))',modifications)
            #get all modifications

        modpos = re.findall(r'[A-Z]([0-9]+)',modifications)
        # get all modification positions
        modpos = [int(d)-1 for d in modpos] #convert all to integers
        if 'N-Term' in modifications:
            modpos.insert(0,0)
            # first of the modkeys will be the N-Term modification
        modi_len = len(modpos)                
        mod_dict = defaultdict(list)
        for (key,value) in zip(modpos, modkeys):
            mod_dict[key].append(value)
        for key in mod_dict:
            mod_dict[key] = sorted(mod_dict[key])

        for ix, s in enumerate(sequence.upper()):
            seqmodi += s
            #print(ix, seqmodi)
            if ix in mod_dict:
                if s == 'X':  # deal with the X first
                    to_replace = [x for x in mod_dict[ix] if x in amino_acids]
                    
                    if len(to_replace) == 1:
                        #print(seqlist[ix],to_replace)
                        if seqlist[ix].islower():
                            seqlist[ix] = to_replace[0][1].lower()
                        else:
                            seqlist[ix] = to_replace[0][1]
                        modi_len += -1  # don't count X == amino acid as a
                        #modification since technically it is not a PTM
                        #mod_dict[ix].remove(to_replace)  # not sure if we want
                        #to remove this or not, but we certainly can if we want
                    elif len(to_replace) == 0:
                         pass  # not an amino acid (listed above at least)
                    #Probably an unidentified mass of the form X10(110.0)
                    
                    else:
                        print('Error parsing sequence {} with '\
                              'modifications {}'.format(sequence, modifications))
                for modi in mod_dict[ix]:
                    if modi in modtext:
                        seqmodi += modtext[modi]
                    elif modi in amino_acids:
                        seqmodi += modi  # for the X amino acids
                    else:
                        logging.warning('New modification {} that'\
                                        ' is not found in sequence {}'.format(
                                             modi, sequence))                

    sequence = ''.join(seqlist)
    return sequence, seqmodi, modi_len, label
    
def pept_print(df,usrdata):
     # Here, try this - returns peptides that are present within the data
    '''
    Lookup info from userdata. Also sorts based on alphabet. 
    '''

    IDquery = df['gene_GeneID']
    #print(IDquery)  # for debugging
    try : 
        matches = usrdata[usrdata.psm_GeneID == IDquery]
        # need to do check for if matches is empty
        #print('matches : {}'.format(matches))

        matches.Sequence = matches.Sequence.astype('object').str.upper()
        #matches.drop_duplicates(cols='Sequence', take_last=False,
        #inplace = True) # cols is depreciated, use subset    
        matches.drop_duplicates(subset='Sequence', keep='first',
                                inplace=True)

        protcount = matches.Sequence.count() # counting peptides (not proteins)
        protcount_S = matches[matches['psm_PSM_IDG' ] < 4].Sequence.count()
        uniques = matches[matches.psm_GeneCount == 1].Sequence.count()
        uniques_S = matches[(matches.psm_GeneCount == 1) &
                            (matches['psm_PSM_IDG' ] < 4)].Sequence.count()
        pepts_str = '_'.join(sorted(matches.Sequence))
    except AttributeError as e: #have never had this happen. Shouldn't occur
         #since the df is derived from usrdata in the first place
        print(e)        
        pepts_str = ''
        protcount = 0
        uniques = 0
    return (set(sorted(matches.Sequence)), pepts_str, protcount, uniques,
         protcount_S, uniques_S)        
    
def e2g_PSM_helper(gene_df_ID, data,EXPTechRepNo):
    total = sum(data[
                data['psm_GeneID']==gene_df_ID]
                ['psm_PSM_useflag'])/EXPTechRepNo
    total_u2g = sum(data[(data['psm_GeneID']==gene_df_ID) &
                         (data['psm_GeneCount' ] == 1)
                        ]['psm_PSM_useflag'])/EXPTechRepNo

    total_S = sum(data[(data['psm_GeneID']==gene_df_ID) &
                       (data['psm_PSM_IDG' ] < 4)
                  ]['psm_PSM_useflag'])/EXPTechRepNo
    total_S_u2g = sum(data[(data['psm_GeneID']==gene_df_ID) &
                           (data['psm_PSM_IDG' ] < 4) &
                           (data['psm_GeneCount' ] == 1)
                      ]['psm_PSM_useflag'])/EXPTechRepNo
    #for match in matches:
    return total, total_u2g, total_S, total_S_u2g
    
   
def area_calculator(gene_df, usrdata, EXPTechRepNo,area_col, normalization):

    matches  = usrdata[(usrdata['psm_GeneID'] == gene_df['gene_GeneID']) &
                           (usrdata['psm_AUC_useflag']==1)] [
                                [area_col, 'psm_GeneCount',
                                 'psm_PSM_IDG', 'MissedCleavages']]
 

    #else : 
       # return None, None, None, None, None, None
        #print('{} : Error - EXPQuantSource is not defined correctly.'.format(
             #datetime.now()))
        #sys.exit(1)
        
    uniq_matches = matches[matches['psm_GeneCount']==1]
    uniq_matches_0 = uniq_matches[uniq_matches['MissedCleavages']==0]
    matches_strict = matches[matches['psm_PSM_IDG'] < 4] 


    values_max = nan_popper([value for value in
                             matches[area_col].values])

    values_adj = nan_popper([value/count for value,count in
                             matches[[area_col,'psm_GeneCount']
                             ].values])
    uniq_values_adj = nan_popper([value/count for value,count in
                                  uniq_matches[[area_col,
                                                'psm_GeneCount']].values])
    uniq_values_adj_0 = nan_popper([value/count for value,count in
                                    uniq_matches_0[[area_col,
                                                    'psm_GeneCount']].values])
    result = (sum(values_max)/normalization,  sum(values_adj)/normalization,
              sum(uniq_values_adj_0)/normalization,
              sum(uniq_values_adj)/normalization)
    return result

def AUC_distributor(inputdata,genes_df,area_col):
     
    #if EXPQuantSource == 'AUC':
    #    inputvalue = inputdata['Precursor Area'] 
    #elif EXPQuantSource == 'Intensity':
    #    inputvalue = inputdata['Intensity']
    #else:
    #    print('{} : Error - EXPQuantSource is not defined correctly.'.format(
    #         datetime.now()))
    #    sys.exit(1)
    inputvalue = inputdata[area_col]
    u2gPept = genes_df[genes_df['gene_GeneID']==inputdata['psm_GeneID']
    ]['gene_GeneArea_gpcAdj_u2g_all'].values
    
    if len(u2gPept) == 1: u2gPept = u2gPept[0] # grab u2g info, should always be
    #of length 1
    elif len(u2gPept) > 1 : 
        print('{} Error - distArea is not singular at GeneID : {}'.format(
             datetime.now(),inputdata['psm_GeneID']))
        # this should never happen (and never has)
    else : 
        distArea = 0
        print('No distArea for GeneID : {}'.format(inputdata['psm_GeneID']))  
    if u2gPept != 0 : 
        totArea = 0
        gene_list = inputdata.psm_GeneList.split(',')
        totArea = sum(genes_df[
             genes_df['gene_GeneID'].isin(gene_list)
                              ].gene_GeneArea_gpcAdj_u2g_all)
        distArea = (u2gPept/totArea) * inputvalue
        #ratio of u2g peptides over total area

    elif u2gPept == 0:  # no uniques, normalize by genecount
        try: distArea = inputvalue/inputdata.psm_GeneCount   
        except ZeroDivisionError: distArea = inputvalue

    return distArea
    
def gene_AUC_sum(genes_df,temp_df,normalization):
    gene = genes_df['gene_GeneID']
    if genes_df['gene_IDSet'] in [1,2]:
        return sum(temp_df[temp_df[
             'psm_GeneID'
        ] == gene].psm_PrecursorArea_dstrAdj)/normalization
    elif genes_df['gene_IDSet'] == 3:
        return 0

def gene_setter(genes_df,genes_df_all, temp_df ):
    IDquery, Pquery = genes_df[['gene_GeneID','gene_PeptideSet']]
    
       # Pquery = set(Pquery.split('_'))
    if any(temp_df[temp_df.psm_GeneID==IDquery]['psm_GeneCount']==1):
       # genes_df.at[pos,'GeneSet'] = 1
        setValue = 1
    elif not any(genes_df_all.gene_PeptideSet.values > Pquery):
        #genes_df.at[pos,'GeneSet'] = 2
        setValue = 2
        #genes_df.ix[genes_df.GeneID == IDquery,'GeneSet'] = 2
    elif any(genes_df_all.gene_PeptideSet.values > Pquery):
        setValue = 3
       # genes_df.at[pos,'GeneSet'] = 3
    try: GeneIDGroup = min(temp_df[temp_df.psm_GeneID == IDquery].psm_PSM_IDG)
    except ValueError: GeneIDGroup = 0
    try : GeneIDGroup_u2g  = min(temp_df[(temp_df.psm_GeneID == IDquery) & (temp_df.psm_GeneCount == 1)].psm_PSM_IDG)
    except ValueError : GeneIDGroup_u2g = 0   
       
    return setValue, GeneIDGroup, GeneIDGroup_u2g   

greaterthan = lambda x, last : True if x < last else False        
def GPG_helper(idset,peptideset, df_all,last):
    # do stuff
    #indx = numgen.next()
    
    #print idset
    if idset == 3:
        gpg =  ''
    else:        
 
        if idset == 1:    
            gpg =  last + 1
        elif idset == 2:
            selection = df_all[df_all.gene_PeptideSet.values == peptideset].gene_GPGroup
            selection = selection[selection!='']
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
                        print('Warning, more than one IDSet type 2 genes '\
                              'already has assignment')
                        gpg =  genes_df['gene_GPGroup'] 
                     
    return gpg, greaterthan(gpg, last)    

    

def GPG_all_helper(genes_df,df_all):

    GPGall = set()
    for pept in genes_df.gene_PeptideSet:
        shared_values = [value for value in
                         df_all[df_all['gene_PeptidePrint'].str.contains(
                              pept,regex=False)].gene_GPGroup.values]
        # regex=False since we don't need regex, and is faster.
        for s in shared_values: 
            if s !='' : GPGall.add(s)
        
    return str([k for k in GPGall]).strip('[').strip(']')    
        #df_all[df_all['_e2g_PeptideSet'].str.contains(pept)
        
def capacity_grabber(geneid, gene_metadata):
    genefraglengths=[]

    for data in gene_metadata[geneid]:
         genefraglengths.append(data[3])
    #sel = df[df.psm_GeneID == geneid][['psm_tProteinList',
    #                                     'psm_ProteinCapacity']]
    #capacity = 0    
    #if len(sel) > 0:
    #    sel.reset_index(inplace = True)
    #    lst_indx = sel.loc[0]['psm_tProteinList'].split(',').index(geneid)
    #    capacity = sel.loc[0]['psm_ProteinCapacity'].split(',')[lst_indx]
    if not genefraglengths:
         genefraglengths.append(0)
         
    return mean(genefraglengths)
        
  
def regex_pattern_all(items):
     patterns=[]
     for item in items:
          pattern = ''
          for element in item:
               pattern += '(?=.*{})'.format(element)
          patterns.append(pattern)

     return patterns

def gene_to_taxon(gene, d):
     if gene:
          gene = str(int(gene))  # just in case
     return d.get(gene)
