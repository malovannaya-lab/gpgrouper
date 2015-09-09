import re, sys
import os
import csv
import numpy as np
from datetime import datetime
from collections import namedtuple, defaultdict
from configparser import SafeConfigParser
from statistics import mean
import logging
from collections import defaultdict
try : 
     from PIL import Image, ImageFont, ImageDraw
    #imagetitle = True
except ImportError: pass
    #imagetitle = False
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
    
    parser = SafeConfigParser()
    cancel = False
    if os.path.isfile('py_config.ini'):
        while True:
            del_confirm = input('A config file is already detected? Are you sure you want to delete and re-run the setup(Y/n)? ')
            if del_confirm == 'Y' :
                os.remove('py_config.ini')
                break
            elif 'n' in del_confirm.lower() : cancel = True
    if not cancel:    
        print('Welcome to PyGrouper setup wizard. First we need to get a list of all taxons you wish to be able to search for.')
        print('Note this requires that you have refseq or peptidome data for each taxon you list.')
        taxons = set()
        while True:
            try:
                if len(taxons) > 0 : print('\nCurrently added taxons : {}\n'.format(' '.join(taxons)))
                newtaxon = input('Enter a taxon ID to add (human is 9606) or Ctrl+C if done : ')
                taxons.add(newtaxon)
            except KeyboardInterrupt:
                break   
        parser.add_section('refseq locations')
        parser.add_section('refseq file sizes')        
        for taxon in taxons:
            while True:
                taxon_loc = input('\nEnter location of the refseq file for taxon {} : '.format(taxon))
                if os.path.isfile(taxon_loc):
                    print('Setting {} refseq location to {}'.format(taxon, taxon_loc))
                    parser.set('refseq locations',taxon,taxon_loc)
                    print('Calculating...')
                    parser.set('refseq file sizes',taxon, str(bufcount(taxon_loc)))
                    print('Done!\n')
                    break
                else : 
                    print('\nFile not found!')
                    if input("Type 'pass' to skip taxon {} from list or [Enter] to try again : ".format(taxon)).lower() == 'pass':
                        break
        parser.write(open('py_config.ini','w'))
        print('Configuration of PyGrouper is complete.\n')
    elif cancel: print('Canceling..\n')
   
def simple_generator_function(size):

    for i in range(size): 
        yield i

def enhanced_generator_counter():
    i = 0
    while True:
        val = (yield i)
        # If value provided, change counter
        if val is not None:
            i = val + 1
        else:
            i += 1    

def csv_reader(inputfile,sep='\t'):
    FastaRecord = namedtuple('FastaRecord','taxonid, homologeneid, geneid, proteingi, fasta')
    with open(inputfile, mode ='rU') as f:
        reader = csv.DictReader(f,delimiter = sep)
        for line in reader:
            record = FastaRecord._make([line['faa_TaxonID'], line['faa_HomologeneID'], line['faa_GeneID'], line['faa_ProteinGI'], line['faa_FASTA']])
            yield record


def file_reader(inputfile,sep='\t'):
     with open(inputfile, mode ='rU') as f:
        for line in f:
            yield namedtuple('taxonid')
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
            elif len(seq) == 0 and chop : frags.append(chop) #special case for if last amino acid is one of the cut sites
            break #no more cuts, wrap up and go home
        if seq[0] not in exceptions:
            frags.append(chop)
            #print(chop, seq)
            chop = ''
        
    frags_1, frags_2, frags_3 = [], [], []
    for i in range(len(frags)-1):
        frags_1.append(frags[i]+frags[i+1])
    for i in range(len(frags)-2):
        frags_2.append(frags[i]+frags[i+1]+frags[i+2])
    for i in range(len(frags)-3):
        frags_2.append(frags[i]+frags[i+1]+frags[i+2]+frags[i+3])    
    #print frags, frags_1, frags_2
    [a.append(a[0][1:]) for a in [frags, frags_1, frags_2, frags_3] if len(a) > 0] # chop off each N-term
    frags = [x for x in frags if len(x) >= minlen]
    fragments = [x for x in frags + frags_1 + frags_2 if len(x) >= minlen]  
    return fragments, len(frags) #len(frags) is no miscuts
    
def protease_old(seq,minlen = 0,cutsites = [], exceptions = []):
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
        if cut == 0 or len(seq) == 0: 
            if len(seq) != 0 :frags.append(seq)
            break #no more cuts, wrap up and go home
        if seq[0] not in exceptions:
            frags.append(chop)
            chop = ''
        #print chop, seq
    frags_1, frags_2, frags_3 = [], [], []
    for i in range(len(frags)-1):
        frags_1.append(frags[i]+frags[i+1])
    for i in range(len(frags)-2):
        frags_2.append(frags[i]+frags[i+1]+frags[i+2])
    for i in range(len(frags)-3):
        frags_2.append(frags[i]+frags[i+1]+frags[i+2]+frags[i+3])    
    #print frags, frags_1, frags_2
    [a.append(a[0][1:]) for a in [frags, frags_1, frags_2, frags_3] if len(a) > 0] # chop off each N-term
    frags = [x for x in frags if len(x) >= minlen]
    fragments = [x for x in frags + frags_1 + frags_2 if len(x) >= minlen]  
    return fragments, len(frags) #len(frags) is no miscuts

def trypsin(seq,minlen):
    frags = []
    chop = ''
    while seq:
        k, r = seq.find('K'), seq.find('R')
        cut = min(k, r)+1 if k >= 0 and r >= 0 else max(k, r)+1
        chop += seq[:cut]
        seq = seq[cut:] 
        #print chop, seq
        if cut == 0 or len(seq) == 0: 
            if len(seq) != 0 :frags.append(seq)
            break #no more cuts, wrap up and go home
        if seq[0] != 'p':
            frags.append(chop)
            chop = ''
    frags_1, frags_2, frags_3 = [], [], []
    for i in range(len(frags)-1):
        frags_1.append(frags[i]+frags[i+1])
    for i in range(len(frags)-2):
        frags_2.append(frags[i]+frags[i+1]+frags[i+2])
    for i in range(len(frags)-3):
        frags_2.append(frags[i]+frags[i+1]+frags[i+2]+frags[i+3])    
    #print frags, frags_1, frags_2
    [a.append(a[0][1:]) for a in [frags, frags_1, frags_2, frags_3] if len(a) > 0] # chop off each N-term
    return [x for x in frags + frags_1 + frags_2 if len(x) >= minlen]     
    
def genematcher(seq, metadata, prot):
    seq = seq.upper()
    #print(type(metadata))
    if not isinstance(metadata,list):
        metadata = []
    if seq in prot:
        for data in prot[seq]:
            metadata.append(data)
    return metadata

def meta_extractor(metadata):
    #taxonid = set()
    genes = defaultdict(list)
    homologenes = set()
    proteingis = set()
    genefraglen = []
    
    if metadata:
        for data in metadata:
            genes[data.geneid].append((data.taxonid, data.homologeneid, data.proteingi, data.genefraglen))
        for gene in genes:
            frags = []
            for genedata in genes[gene]:
                #taxonid.add(genedata[0])
                homologenes.add(genedata[1])
                proteingis.add(genedata[2])
                frags.append(genedata[3])
            if len(frags) > 0 : genefraglen.append(str(mean(frags)))
            else : genefraglen.append(0)
    genelst = [gene for gene in genes]
    return ','.join(genelst), ','.join(proteingis), len(genelst), len(proteingis), ','.join(homologenes), len(homologenes), ','.join(genefraglen)             


def genematcher_leg(seq,genelist, proteinlist, genecount, proteincount, metadata, prot):
    ''' matches genes via refseq in prot dictionary, which has a list of these guys : namedtuple('RefseqInfo','geneid, homologeneid, proteingi,genefraglen')
    '''
    seq = seq.upper()
    # for key in prot.keys():
        # if seq in key:
    if genecount == 0:    
        genelist = []
        proteinlist  = []
        metadata = []
        #homologeneid = []
        #proteincapacity = defaultdict(list)
        
    else:
        genelist = genelist.split(',')
        proteinlist = proteinlist.split(',')
        #homologeneid = homologeneid.split(',')
        #proteincapacity = proteincapacity.split(',')
    if seq in prot:    
        for lst in prot[seq]:
            genelist.append(lst.geneid)
            proteinlist.append(lst.proteingi)
            metadata

            homologeneid.append(lst.homologeneid)
            proteincapacity.append(str(lst.genefraglen))
    #except KeyError:
    #    pass #might happen for some peptides. A log of these guys is noted.
    #except AttributeError as e:
    #    print(e) #have never seen this happen
    return ','.join(genelist), ','.join(proteinlist), len(genelist), len(proteinlist), ','.join(homologeneid), ','.join(proteincapacity)

        
def dict_modifier(d,exceptions = {},exception_float = True):
    """Dictionary modifier function for the user. Takes in dictionary (d), as well as all dictionary keys that do not have to be floats (exceptions),
    the valid strings that the exception key values can hold (valid_exceptions), as well as a boolean value for if the exception key values can be numbers
    """
    while True:
        print('\nThe current settings are : ')
        for key,value in list(d.items()): print('{} : {}'.format(key,value))
        try:
            
            filter_change = input('Type a value to change (case sensitive), or hit Ctrl+C if done. ')
            try :
                new_value = input('The current value is set to {}. Type a new value. : '.format(d[filter_change]))
                if filter_change not in list(exceptions.keys()):
                    try : 
                        new_value_float = int(new_value)
                        logging.info('{} changed from {} to {}'.format(filter_change, d[filter_change], new_value_float))
                        d[filter_change] = new_value_float   
                    except ValueError:
                        print('\nInvalid entry. New value must be an integer.\n')
                elif filter_change in list(exceptions.keys()):
                    if new_value in exceptions[filter_change]:
                        logging.info('{} changed from {} to {}'.format(filter_change, d[filter_change], new_value))
                        d[filter_change] = new_value
                        
                    elif exception_float:
                        try : 
                            new_value_float = float(new_value)
                            logging.info('{} changed from {} to {}'.format(filter_change, d[filter_change], new_value_float))
                            d[filter_change] = float(new_value_float)
                        except ValueError:
                            print('\nInvalid entry. New value must be a number or in {}.\n'.format(valid_exceptions))
                    else : print('\nInvalid entry. New value must be in {}.\n'.format(valid_exceptions))
            except KeyError :
                print('\nInvalid filter value.\n')
                
        except KeyboardInterrupt:
            print('')
            break
    return d   

def nanPopper(lst):
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
    if df['Charge'] < d['Filter_Z_min'] or df['Charge'] > d['Filter_Z_max'] : # if charge is empty (nan), this will not be true
        AUC_flag = 0
        PSM_flag = 0

    elif df['_data_SequenceModiCount'] > d['Filter_Modi'] : 
        AUC_flag = 0
        PSM_flag = 0

    elif np.isnan(df['IonScore']) and np.isnan(df['q-Value']): 
        AUC_flag = 1 #retains WLs Q2up assignments
        PSM_flag = 0

    elif df['IonScore'] < d['Filter_IS'] : 
        AUC_flag = 0
        PSM_flag = 0
    elif df['q-Value'] > d['Filter_qV'] : 
        AUC_flag = 0
        PSM_flag = 0
    elif df['PEP'] > d['Filter_PEP'] : 
        AUC_flag = 0
        PSM_flag = 0
    elif df['_data_PSM_IDG'] > d['Filter_IDG'] : 
        AUC_flag = 0
        PSM_flag = 0
    elif df['_data_PeptRank'] == 0 : 
        AUC_flag = 0 # omit multiPSM peak areas
        PSM_flag = 1  # Change from 0 in grouper
    else: 
        AUC_flag = 1
        PSM_flag = 1
    return AUC_flag, PSM_flag
    
def PSM_flagger(df):
    if df['IonScore']=='' and df['q-Value'] =='': PSM_flag = 0 # omit WL's Q2up peptides
    #elif df['data_AUC_nUseFLAG'] ==0: PSM_flag = 0 #omit not used for area, note not the same as no area good PSMs
    elif Filter_IS>0 and df['IonScore'] < Filter_IS : PSM_flag = 0
    elif Filter_qV > 0 and df['q-Value'] > Filter_qV : PSM_flag = 0
    elif Filter_PEP > 0 and df['PEP'] > Filter_PEP : PSM_flag = 0
    elif Filter_IDG > 0 and df['data_PSM_IDG'] > Filter_IDG : PSM_flag= 0
    elif Filter_Z_min > 0 and df['Charge'] < Filter_Z_min : PSM_flag = 0
    elif Filter_Z_max > 0 and df['Charge'] > Filter_Z_max : PSM_flag = 0
    elif Filter_Modi > 0 and df['_data_SequenceModiCount'] > Filter_Modi : PSM_flag= 0
    else: PSM_flag = 1
    return PSM_flag   
#list(usrdata.loc[k, ['IonScore','q-Value','data_PeptRank','Charge','SeqModiCount','PEP','data_PSM_IDG']])+[Filter_qV,Filter_IS,Filter_PEP,Filter_IDG,Filter_Z_min,Filter_Z_max,Filter_Modi]

def IDG_picker_old(df):
   # try: 

    if (df['IonScore']>= 30 and df['q-Value'] <= .01): IDGout = 1
    elif (df['IonScore'] >= 30 and df['q-Value'] <= .05): IDGout = 2
    elif (df['IonScore'] >= 20 and df['q-Value'] <= .01): IDGout = 3
    elif (df['IonScore'] >= 20 and df['q-Value'] <= .05): IDGout = 4
    elif (df['IonScore'] >= 10 and df['q-Value'] <= .01): IDGout = 5       
    elif (df['IonScore'] >= 10 and df['q-Value'] <= .05): IDGout = 6
    elif (df['IonScore'] >= 0 and df['q-Value'] <= .01): IDGout = 7
    elif (df['IonScore'] >= 0 and df['q-Value'] <= .05): IDGout = 8
    else: IDGout  = 9
   # except: IDGout = 9 #special case for when ion score or qvalue is missing
    return IDGout
 
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
 
def SeqModi_old(df):
    modtext = {'(DeStreak)' : '(des)', '(Deamidated)' : '(dam)', '(Carbamidomethyl)' : '(car)', '(Oxidation)' : '(oxi)', '(Phospho)' : '(pho)',
    '(Acetyl)': '(ace)', '(GlyGly)' : '(gg)', '(Label:13C(6))' : '(lab)', '(Label:13C(6)+GlyGly)' : '(labgg)'}
        
    if any(c for c in df['Sequence'] if c.islower()):
        if 'Label' in df['Modifications'] : label = 1 # for indication of label for peptide
        else : label = 0
        modpos = re.findall(r'([0-9]+)',df['Modifications'])
        modpos = [int(d) for d in modpos]
        modkeys = re.findall(r'(\(\w+\))',df['Modifications'])
        if 'N-Term' in df['Modifications']: modpos.insert(0,1)
        querylst = list(df['Sequence'])
        for k in range(len(modpos)):
            try:
                if k != 0:
                    if modpos[k] != modpos[k-1]: querylst[modpos[k]-1] = querylst[modpos[k]-1].upper()+modtext[modkeys[0]]
                    else : querylst[modpos[k]-1] = querylst[modpos[k]-1]+modtext[modkeys[0]]
                else :  querylst[modpos[k]-1] = querylst[modpos[k]-1].upper()+modtext[modkeys[0]]
            except Exception as e:  #due to key error  
                #print e
                querylst[modpos[k]-1] = modkeys[0]
                    
              
            modkeys.pop(0) # pop the first guy and go again

        return ''.join(querylst), len(modpos), label
    else: 
        label = 0 # for indication of label of peoptide
        return df['Sequence'], 0, label 

def SeqModi_nope(sequence, modifications):
    modtext = {'(DeStreak)' : '(des)', '(Deamidated)' : '(dam)', '(Carbamidomethyl)' : '(car)', '(Oxidation)' : '(oxi)', '(Phospho)' : '(pho)',
    '(Acetyl)': '(ace)', '(GlyGly)' : '(gg)', '(Label:13C(6))' : '(lab)', '(Label:13C(6)+GlyGly)' : '(labgg)'}
        
    if any(c for c in sequence if c.islower()):
        if 'Label' in modifications : #routine for grabbing everything enclosed in parenthesis, including nested. Easier to write than regex when dealing with nested
            label = 1 # for indication of label for peptide
            openb, modi, modkeys = 0, '', []
            for e in modifications:
                if e == '(':
                    openb += 1
                elif e == ')':
                    openb += -1
                if openb:
                    modi += e
                elif not openb:
                    if modi:
                        modi += e #add the last parenthesis
                        modkeys.append(modi)
                    modi = ''
                
        else :
            label = 0
            modkeys = re.findall(r'(\([^\)]+\))',modifications)
        modpos = re.findall(r'[A-Z]([0-9]+)',modifications) #only numbers after letters
        modpos = [int(d) for d in modpos]
        
        if 'N-Term' in modifications: modpos.insert(0,0)
        querylst = list(sequence)
        for k in range(len(modpos)):
            try:
                if k != 0:
                    if modpos[k] != modpos[k-1]: 
                        querylst[modpos[k]-1] = querylst[modpos[k]-1].upper()+modtext[modkeys[0]]
                        #print(''.join(querylst))
                    else : 
                        querylst[modpos[k]-1] = querylst[modpos[k]-1]+modtext[modkeys[0]]
                        #print(''.join(querylst))
                else :  
                    querylst[modpos[k]] = querylst[modpos[k]].upper()+modtext[modkeys[0]]
                    #print(''.join(querylst))
            except KeyError:  # skip modification
                #print('KeyError : {}'.format(modkeys[0]))
                pass
               # querylst[modpos[k]-1] = modkeys[0]
            except Exception as e:
                print(e)
                print(sequence, modifications)
              
            if modkeys : modkeys.pop(0) # pop the first guy and go again

        return ''.join(querylst), len(modpos), label
    else: 
        label = 0 # for indication of label of peoptide
        return sequence, 0, label  
        
def seq_modi(sequence, modifications):

    '''
    function to output a modified string of sequence that includes all modifications at the appropriate place.
    Modified amino acids are indicated in the input 'sequence' by a lowercase letter.
    A dictionary of potential modification keys and their corresponding substitutions in the primary sequence is provided;
    all other modifications should be ignored.
    Sometimes N-terminal modifications are present, this should be considered the same as a modification at the first amino acid.
    '''
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    amino_acids = ['('+x+')' for x in amino_acids]
    modtext = {'(DeStreak)' : '(des)', '(Deamidated)' : '(dam)', '(Carbamidomethyl)' : '(car)', '(Oxidation)' : '(oxi)', '(Phospho)' : '(pho)',
    '(Acetyl)': '(ace)', '(GlyGly)' : '(gg)', '(Label:13C(6))' : '(lab)', '(Label:13C(6)+GlyGly)' : '(labgg)'}
    seqmodi = ''
    seqlist = list(sequence)
    modi_len = 0  # default length is zero, can change
    label = 0  # default label = 0, can change
    if any(c for c in sequence if (c.islower() or c == 'X')): # check if any modifications to deal with
        if 'Label' in modifications: #routine for grabbing everything enclosed in parenthesis, including nested. Easier to write than regex when dealing with nested
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
            modkeys = re.findall(r'(\([^\)]+\))',modifications)  #get all modifications

        modpos = re.findall(r'[A-Z]([0-9]+)',modifications)  # get all modification positions
        modpos = [int(d)-1 for d in modpos] #convert all to integers
        if 'N-Term' in modifications:
            modpos.insert(0,0)  # first of the modkeys will be the N-Term modification
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
                        modi_len += -1  # don't count X == amino acid as a modification since technically it is not a PTM
                        #mod_dict[ix].remove(to_replace)  # not sure if we want to remove this or not, but we certainly can if we want
                    elif len(to_replace) == 0:
                         pass  # not an amino acid (listed above at least) Probably an unidentified mass of the form X10(110.0)
                    
                    else:
                        print('Error parsing sequence {} with modifications {}'.format(sequence, modifications))
                for modi in mod_dict[ix]:
                    if modi in modtext:
                        seqmodi += modtext[modi]
                    elif modi in amino_acids:
                        seqmodi += modi  # for the X amino acids
                    else:
                        logging.warning('New modification {} that is not found in sequence {}'.format(modi, sequence))                

    sequence = ''.join(seqlist)
    return sequence, seqmodi, modi_len, label
    
    
def seq_modi_old(sequence, modifications):

    '''
    function to output a modified string of sequence that includes all modifications at the appropriate place.
    Modified amino acids are indicated in the input 'sequence' by a lowercase letter.
    A dictionary of potential modification keys and their corresponding substitutions in the primary sequence is provided;
    all other modifications should be ignored.
    Sometimes N-terminal modifications are present, this should be considered the same as a modification at the first amino acid.
    '''

    modtext = {'(DeStreak)' : '(des)', '(Deamidated)' : '(dam)', '(Carbamidomethyl)' : '(car)', '(Oxidation)' : '(oxi)', '(Phospho)' : '(pho)',
    '(Acetyl)': '(ace)', '(GlyGly)' : '(gg)', '(Label:13C(6))' : '(lab)', '(Label:13C(6)+GlyGly)' : '(labgg)'}
    seqmodi = ''
    seqlist = list(sequence)
    modi_len = 0  # default length is zero, can change
    label = 0  # default label = 0, can change
    if any(c for c in sequence if (c.islower() or c == 'X')): # check if any modifications to deal with
        if 'Label' in modifications: #routine for grabbing everything enclosed in parenthesis, including nested. Easier to write than regex when dealing with nested
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
            
            modkeys = re.findall(r'(\([^\)]+\))',modifications)  #get all modifications

        modpos = re.findall(r'[A-Z]([0-9]+)',modifications)  # get all modification positions
        modpos = [int(d)-1 for d in modpos] #convert all to integers

        if 'N-Term' in modifications:
            if modpos: 
                 if modpos[0] == 0:
                      modkeys[0] = modkeys[0] + '_' + modkeys[1]
                      modkeys.pop(1)
                 else:
                      modpos.insert(0,0)
            else:
                modpos.insert(0,0)
        if any(modpos.count(x) > 1 for x in modpos):
            for ix, x in enumerate(modpos):
                if modpos.count(x) > 1 and modpos.index(x) == ix: #just once per x
                    dups =  [ix for ix, y in enumerate(modpos) if y == x]
                    modkeys[dups[0]] = '_'.join(sorted(modkeys[z] for z in dups))
                    pops = [z for z in dups if z != dups[0]]
                    for p in pops: modkeys.pop(p)

        #print(modkeys, modpos)
        modi_len = len(modpos)
        print(modkeys)
        for ix, s in enumerate(sequence.upper()):
            print(seqmodi)
            seqmodi += s
            #print(ix, seqmodi)
            if ix in modpos:
                if modkeys[modpos.index(ix)] in modtext:
                    seqmodi += modtext[modkeys[modpos.index(ix)]]

                elif modkeys[modpos.index(ix)] not in modtext:
                    if s == 'X':
                        seqmodi += modkeys[modpos.index(ix)]
                        seqlist[ix] = modkeys[modpos.index(ix)][1] #replace sequence with the predicted amino acid
                    else:
                        mods = modkeys[modpos.index(ix)].split('_') # deal with multiple modifications on single amino acid
                        for mod in mods:
                            if mod in modtext:
                                seqmodi += modtext[mod]
                            else:
                                print('New modification that is not found : {}'.format(modkeys[modpos.index(ix)]))
                                logging.warning('New modification that is not found : {}'.format(modkeys[modpos.index(ix)]))


    sequence = ''.join(seqlist)
    return sequence, seqmodi, modi_len, label

def SeqModi_old123(sequence, modifications):
    modtext = {'(DeStreak)' : '(des)', '(Deamidated)' : '(dam)', '(Carbamidomethyl)' : '(car)', '(Oxidation)' : '(oxi)', '(Phospho)' : '(pho)',
    '(Acetyl)': '(ace)', '(GlyGly)' : '(gg)', '(Label:13C(6))' : '(lab)', '(Label:13C(6)+GlyGly)' : '(labgg)'}
        
    if any(c for c in sequence if c.islower()):
        if 'Label' in modifications : #routine for grabbing everything enclosed in parenthesis, including nested. Easier to write than regex when dealing with nested
            label = 1  # for indication of label for peptide
            openb, modi, modkeys = 0, '', []
            for e in modifications:
                if e == '(':
                    openb += 1
                elif e == ')':
                    openb += -1
                if openb:
                    modi += e
                elif not openb:
                    if modi:
                        modi += e #add the last parenthesis
                        modkeys.append(modi)
                    modi = ''
                
        else :
            label = 0
            modkeys = re.findall(r'(\([^\)]+\))',modifications)
        modpos = re.findall(r'[A-Z]([0-9]+)',modifications) #only numbers after letters
        modpos = [int(d) for d in modpos]
        
        if 'N-Term' in modifications: modpos.insert(0,0)
        querylst = list(sequence)
        for k in range(len(modpos)):
            try:
                if modpos[k] != 0:
                    if modpos[k] != modpos[k-1]: querylst[modpos[k]-1] = querylst[modpos[k]-1].upper()+modtext[modkeys[0]]
                    else : querylst[modpos[k]-1] = querylst[modpos[k]-1].upper()+modtext[modkeys[0]]
                else :  querylst[modpos[k]] = querylst[modpos[k]].upper()+modtext[modkeys[0]]
            except KeyError:  # skip modification
                pass
               # querylst[modpos[k]-1] = modkeys[0]
            except IndexError:
                querylst[modpos[k]-1] = querylst[modpos[k]-1].upper()+modtext[modkeys[0]]
            except Exception as e:
                print(e)
                print(sequence, modifications)
            if modkeys : modkeys.pop(0) # pop the first guy and go again

        return ''.join(querylst), len(modpos), label
    else: 
        label = 0 # for indication of label of peoptide
        return sequence, 0, label         
        
        


def GeneFinder1_old(df,peptidome,usrfile):
    '''
    Match the input fragment via index, then extract genelist and protein list from proteome dict.
    '''
    GeneList = set()
    ProteinList = set()
    sequence = df['Sequence'].upper()
    if df._data_GeneCount == 0:
        try:
            genes = peptidome[sequence][1]
            proteins = peptidome[sequence][2]
            for gene in genes: GeneList.add(gene)
            for protein in proteins: ProteinList.add(protein)
        except KeyError:
            pass
            #print 'No match for sequence {} in {}'.format(sequence,usrfile)
           # logging.warning('No match for sequence {} in {}'.format(sequence,usrfile))

        return ','.join(GeneList), ','.join(ProteinList), len(GeneList), len(ProteinList)       
    else:
        return df._data_tGeneList, df._data_tProteinList, df._data_GeneCount, df._data_ProteinCount
 # 0 has the taxonID, 1 has the sequences, 2 has the genelist,    
# 16 has genecount, 17 protein list


def GeneFinder1(sequence, tgenelist, tproteinlist, genecount, proteincount, peptidome):
    '''
    Match the input fragment via index, then extract genelist and protein list from proteome dict.
    '''
    #*args = tgenelist, tproteinlist, genecount, proteincount
    
    
    GeneList = set()
    ProteinList = set()
    sequence = sequence.upper()
    if genecount == 0:
        try:
            genes = peptidome[sequence][1]
            proteins = peptidome[sequence][2]
            for gene in genes: GeneList.add(gene)
            for protein in proteins: ProteinList.add(protein)
        except KeyError:
            pass
            #print 'No match for sequence {} in {}'.format(sequence,usrfile)
           # logging.warning('No match for sequence {} in {}'.format(sequence,usrfile))

        return ','.join(GeneList), ','.join(ProteinList), len(GeneList), len(ProteinList)       
    else:
        return tgenelist, tproteinlist, genecount, proteincount
        #pass
        #return df._data_tGeneList, df._data_tProteinList, df._data_GeneCount, df._data_ProteinCount
 # 0 has the taxonID, 1 has the sequences, 2 has the genelist,    
# 16 has genecount, 17 protein list


def PeptPrint(df,peptidome): # set(sequence.upper()) and sort it. 
    '''
    Lookup info from peptidome table
    '''
# 0 has the taxonID, 1 has the sequences, 2 has the genecount, 3 has the protein list   
    try: 
        pepts = peptidome[df['GeneID']][1]
        protlist = peptidome[df['GeneID']][2]
        pepts_str =  '_'.join(pepts)
        #genecount = len(genelist)
        protcount = len(pepts)
        uniques = 0
        #for number in peptidome[df['GeneID']][2]:
        for number in protlist:
            if number == '1': uniques +=1
        
        #
        #uniques = 0 
       # for pept in pepts: # we don't want to go down this road boys
       #     duplicate = False
       #     for key in peptidome.keys(): 
       #         if key !=df['GeneID']:
       #             if pept in peptidome[key][1] : duplicate = True 
       #             else: uniques += 1
       #         if duplicate : break  
    except KeyError: 
        pepts_str = None
        #genecount = 0
        protcount = 0
        uniques = None
    return pepts_str, protcount, uniques

def PeptPrint1(df,usrdata): # Here, try this - returns peptides that are present within the data
    '''
    Lookup info from userdata. Also sorts based on alphabet. 
    '''
    #try :
    IDquery = df['_e2g_GeneID']
    #print(IDquery)
    try : 
        matches = usrdata[usrdata._data_GeneID == IDquery] # need to do check for if matches is empty
        #print('matches : {}'.format(matches))
     #   print matches
        matches.Sequence = matches.Sequence.astype('object').str.upper()
    #    matches.drop_duplicates(cols='Sequence', take_last=False, inplace = True) # cols is depreciated, use subset    
        matches.drop_duplicates(subset='Sequence', take_last=False, inplace = True)

        protcount = matches.Sequence.count() # counting peptides (not proteins)
        protcount_S = matches[matches['_data_PSM_IDG' ] < 4].Sequence.count()
        
        uniques = matches[matches._data_GeneCount == 1].Sequence.count()
        uniques_S = matches[(matches._data_GeneCount == 1) & (matches['_data_PSM_IDG' ] < 4)].Sequence.count()
        pepts_str = '_'.join(sorted(matches.Sequence))
        # pepts = set()
        # for match in matches.Sequence  : pepts.add(match.upper())  # build a set of peptides (unique)
        # pepts_str = '_'.join(sorted(pepts))  # sorted alphabetically    
        # protcount = len(pepts)
        # uniq_matches = matches[matches.data_GeneCount == 1 ].Sequence
        # uniques = len(uniq_matches) # uniques is just a count
    except AttributeError as e: #have never had this happen. Shouldn't occur since the df is derived from usrdata in the first place
        print(e)        
        pepts_str = ''
        protcount = 0
        uniques = 0
    return set(sorted(matches.Sequence)), pepts_str, protcount, uniques, protcount_S, uniques_S        
    
def e2g_PSM_helper(gene_df_ID, data,EXPTechRepNo):
    total = sum(data[data['_data_GeneID']==gene_df_ID]['_data_PSM_nUseFLAG'])/EXPTechRepNo
    total_u2g = sum(data[(data['_data_GeneID']==gene_df_ID) & (data['_data_GeneCount' ] == 1)  ] ['_data_PSM_nUseFLAG'])/EXPTechRepNo

    total_S = sum(data[(data['_data_GeneID']==gene_df_ID) & (data['_data_PSM_IDG' ] < 4)  ] ['_data_PSM_nUseFLAG'])/EXPTechRepNo
    total_S_u2g = sum(data[(data['_data_GeneID']==gene_df_ID) & (data['_data_PSM_IDG' ] < 4) & (data['_data_GeneCount' ] == 1) ] ['_data_PSM_nUseFLAG'])/EXPTechRepNo
    #for match in matches:
    return total, total_u2g, total_S, total_S_u2g
    
def e2g_PSM_helper_old(gene_df, data,EXPTechRepNo):
    total = sum(data[data['_data_GeneID']==gene_df['_e2g_GeneID']]['_data_PSM_nUseFLAG'])/EXPTechRepNo
    total_u2g = sum(data[(data['_data_GeneID']==gene_df['_e2g_GeneID']) & (data['_data_GeneCount' ] == 1)  ] ['_data_PSM_nUseFLAG'])/EXPTechRepNo
    match  = data[data['_data_GeneID'] == gene_df['_e2g_GeneID']]['_data_GeneID']
    total_S = sum(data[(data['_data_GeneID']==gene_df['_e2g_GeneID']) & (data['_data_PSM_IDG' ] < 4)  ] ['_data_PSM_nUseFLAG'])/EXPTechRepNo
    total_S_u2g = sum(data[(data['_data_GeneID']==gene_df['_e2g_GeneID']) & (data['_data_PSM_IDG' ] < 4) & (data['_data_GeneCount' ] == 1) ] ['_data_PSM_nUseFLAG'])/EXPTechRepNo
    #for match in matches:
    
    # try:
        # if any(value =='1' for value in peptidome[gene_df['GeneID']][2]):
            # total_uniq = sum( data[data['GeneID']==gene_df['GeneID']]['data_PSM_nUseFLAG']) / EXPTechRepNo
    # except KeyError: 
        # print 'KeyError : fail to find GeneID {} within peptidome'.format(gene_df['GeneID'])
        # print gene_df['GeneID']
    return total, total_u2g, total_S, total_S_u2g

   
def AreaCalculator(gene_df, usrdata, EXPTechRepNo,EXPQuantSource):
    if EXPQuantSource == 'AUC':
        matches  = usrdata[(usrdata['_data_GeneID'] == gene_df['_e2g_GeneID']) & (usrdata['_data_AUC_nUseFLAG']==1)] [['Precursor Area','_data_GeneCount','_data_PSM_IDG', '# Missed Cleavages']]
        #matches  = usrdata[(usrdata['GeneID'].isin(gene_df['GeneID'])) & (usrdata['data_AUC_nUseFLAG']==1)] [['Precursor Area','data_GeneCount','data_PSM_IDG', '# Missed Cleavages']]

        normalization = 10**9
 
    elif EXPQuantSource == 'Intensity':
        matches  = usrdata[(usrdata['_data_GeneID'] == gene_df['_e2g_GeneID']) & (usrdata['data_AUC_nUseFLAG']==1)] [['Intensity','_data_GeneCount']]
        normalization = 10**5        

    else : 
       # return None, None, None, None, None, None
        print('{} : Error - EXPQuantSource is not defined correctly.'.format(datetime.now()))
        sys.exit(1)
        
    uniq_matches = matches[matches['_data_GeneCount']==1]
    uniq_matches_0 = uniq_matches[uniq_matches['# Missed Cleavages']==0]
    matches_strict = matches[matches['_data_PSM_IDG'] < 4] 
    #uniq_matches_strict = uniq_matches[usrdata['_data_PSM_IDG'] < 4] 

    values_max = nanPopper([value for value in matches['Precursor Area'].values])
    #values_max_u2g = nanPopper([value for value in matches[matches['_data_GeneCount']==1]['Precursor Area'].values])
    values_adj = nanPopper([value/count for value,count in matches[['Precursor Area','_data_GeneCount']].values])
    uniq_values_adj = nanPopper([value/count for value,count in uniq_matches[['Precursor Area','_data_GeneCount']].values])
    uniq_values_adj_0 = nanPopper([value/count for value,count in uniq_matches_0[['Precursor Area','_data_GeneCount']].values])
    #values_S = nanPopper([value/count for value,count in matches_strict[['Precursor Area','_data_GeneCount']].values])
    #uniq_values_S = nanPopper([value/count for value,count in uniq_matches_strict[['Precursor Area','_data_GeneCount']].values])
    result = sum(values_max)/normalization,  sum(values_adj)/normalization, sum(uniq_values_adj_0)/normalization, sum(uniq_values_adj)/normalization
    return result

def AUC_distributor(inputdata,genes_df,EXPQuantSource):
    if EXPQuantSource == 'AUC':
        inputvalue = inputdata['Precursor Area'] 
    elif EXPQuantSource == 'Intensity':
        inputvalue = inputdata['Intensity']
    else:
        print('{} : Error - EXPQuantSource is not defined correctly.'.format(datetime.now()))
        sys.exit(1)
    u2gPept = genes_df[genes_df['_e2g_GeneID']==inputdata['_data_GeneID']]['_e2g_nGPArea_Sum_u2g_all'].values
    
    if len(u2gPept) == 1: u2gPept = u2gPept[0] # grab u2g info, should always be of length 1
    elif len(u2gPept) > 1 : 
        print('{} Error - distArea is not singular at GeneID : {}'.format(datetime.now(),inputdata['_data_GeneID'])) # this should never happen
    else : 
        distArea = 0
        print('No distArea for GeneID : {}'.format(inputdata['_data_GeneID']))  
    if u2gPept != 0 : 
        totArea = 0
        GeneList = inputdata._data_tGeneList.split(',')
        totArea = sum(genes_df[genes_df['_e2g_GeneID'].isin(GeneList)]._e2g_nGPArea_Sum_u2g_all)
        distArea = (u2gPept/totArea) * inputvalue #ratio of u2g peptides over total area

    elif u2gPept == 0:  # try normalize by genecount
        try:distArea = inputvalue/inputdata._data_GeneCount   
        except ZeroDivisionError: distArea = inputvalue

    return distArea
    
def gene_AUC_sum(genes_df,temp_df,normalization):
    gene = genes_df['_e2g_GeneID']
    if genes_df['_e2g_IDSet'] in [1,2]:
        return sum(temp_df[temp_df['_data_GeneID'] == gene]._data_PrecursorArea_dstrAdj)/normalization
    elif genes_df['_e2g_IDSet'] == 3:
        return 0

def GeneSetter(genes_df,genes_df_all, temp_df ):
    IDquery, Pquery = genes_df[['_e2g_GeneID','_e2g_PeptideSet']]
       # Pquery = set(Pquery.split('_'))
    if any(temp_df[temp_df._data_GeneID==IDquery]['_data_GeneCount']==1):
       # genes_df.at[pos,'GeneSet'] = 1
        setValue = 1
    elif not any(genes_df_all._e2g_PeptideSet > Pquery):
        #genes_df.at[pos,'GeneSet'] = 2
        setValue = 2
        #genes_df.ix[genes_df.GeneID == IDquery,'GeneSet'] = 2
    elif any(genes_df_all._e2g_PeptideSet > Pquery):
        setValue = 3
       # genes_df.at[pos,'GeneSet'] = 3
    try: GeneIDGroup = min(temp_df[temp_df._data_GeneID == IDquery]._data_PSM_IDG)
    except ValueError: GeneIDGroup = 0
    try : GeneIDGroup_u2g  = min(temp_df[(temp_df._data_GeneID == IDquery) & (temp_df._data_GeneCount == 1)]._data_PSM_IDG)
    except ValueError : GeneIDGroup_u2g = 0   
       
    return setValue, GeneIDGroup, GeneIDGroup_u2g
  
  
# def GPGrouper(genes_df,df_all,indx):
    # if genes_df._e2g_IDSet == 3:
        # pass
        
    # else:        
   #     #indx = genes_df.index.values[0] #should be singular
       # print indx
   #     # try :
            # prevG = df_all.ix[indx.next()-1]['_e2g_GPGroup']
            # print prevG
    #    #    prevG, prevG_all = map(lambda x: x, *genes_df.ix[indx][['_e2g_GPGroup','_e2g_GPGroups_All']]) # May return Long value, convert to int at the end
            # prev = True
        # except KeyError:
            # print 'keyerror called'
            # prevG = 0
            # prev = False
     #    #   test, test1 = map(lambda x: x, *df.ix[b][['_e2g_PeptidePrint','_e2g_IDGroup']].values)
        # if genes_df._e2g_IDSet == 1:
            # print prevG
            # print type(prevG)
            # GPGroup = prevG + 1
        # elif genes_df.e2g_IDSet == 2:
            # selection = df_all[df_all._e2g_PeptideSet == genes_df._e2g_PeptideSet]._e2g_GPGroup
            # selection = selection[selection.notnull()]
            # if len(selection) == 1:
                # GPGroup = selection.values[0]
            # elif len(selection) == 0:
                # GPGroup = prevG + 1
            # elif len(selection) > 1:
                # print 'Warning, more than one IDSet type 2 genes already has assignment'
            
            # if not prev : GPGroup = 1 #should never happen due to sort    
        # return GPGroup        
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
            selection = df_all[df_all._e2g_PeptideSet == peptideset]._e2g_GPGroup
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
                        print('Warning, more than one IDSet type 2 genes already has assignment')
                        gpg =  genes_df['_e2g_GPGroup'] 
                     
    return gpg, greaterthan(gpg, last)    

    
def GPG_helper_old(genes_df, df_all,numgen):
    # do stuff
    indx = next(numgen)
    if genes_df['_e2g_IDSet'] == 3:
        return genes_df['_e2g_GPGroup'] 
    else:        
        try :
            prevG = df_all.at[indx-1,'_e2g_GPGroup']
            prev = True
            if prevG == '': prevG = 0
        except KeyError:
            #print 'keyerror called'
            prevG = 0
            prev = False    
        if genes_df['_e2g_IDSet'] == 1:    
            return prevG + 1
        elif genes_df['_e2g_IDSet'] == 2:
            selection = df_all[df_all._e2g_PeptideSet == genes_df['_e2g_PeptideSet']]._e2g_GPGroup
            selection = selection[selection!='']
            if len(selection) == 1:
                return selection.values[0]
            elif len(selection) == 0:
                return prevG + 1
            elif len(selection) > 1:
                uniq_sel = []
                for k in range(len(selection)):
                    if selection.iloc[k] not in uniq_sel: 
                        uniq_sel.append(selection.iloc[k])
                    if len(uniq_sel) == 1:    
                        return uniq_sel[0]
                    elif len(uniq_sel) > 1:
                        print('Warning, more than one IDSet type 2 genes already has assignment')
                        return genes_df['_e2g_GPGroup'] 
                    if not prev : return 1 #should never happen due to sort    
def GPGall_helper(genes_df,df_all):

    GPGall = set()
    for pept in genes_df._e2g_PeptideSet:
        shared_values = [value for value in df_all[df_all['_e2g_PeptidePrint'].str.contains(pept,regex=False)]._e2g_GPGroup.values]
        for s in shared_values: 
            if s !='' : GPGall.add(s)
        
    return str([k for k in GPGall]).strip('[').strip(']')    
        #df_all[df_all['_e2g_PeptideSet'].str.contains(pept)
        
def capacity_grabber(geneid, df):
    sel = df[df._data_GeneID == geneid][['_data_tGeneList','_data_ProteinCapacity']]
    capacity = 0    
    if len(sel) > 0:
        sel.reset_index(inplace = True)
        lst_indx = sel.at[0,'_data_tGeneList'].split(',').index(geneid)
        capacity = sel.at[0,'_data_ProteinCapacity'].split(',')[lst_indx]
    
    return capacity
        
  
