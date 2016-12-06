#===============================================================================#
import re, sys
import os
import csv
import numpy as np
from datetime import datetime
from collections import namedtuple, defaultdict, OrderedDict, deque
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

def byte_formatter(b):
    conv = b/(2**10)
    if conv < 1000:
        return '{:.4f} KB'.format(conv)
    elif conv > 1000:
        conv = conv/(2**10)
        if conv < 1000:
            return '{:.4f} MB'.format(conv)
        elif conv < 1000:
            conv = conv/(2**10)
            return '{:.4f} GB'.format(conv)


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


def rolling_window(seq,length_of_window):
    it = iter(seq)
    window = deque(maxlen=length_of_window)
    for _ in range(length_of_window):
        window.append(next(it))
    yield tuple(window)
    while it:
        try:
            window.append(next(it))
        except StopIteration:
            break
        yield tuple(window)

def protease(seq,minlen = 0, cutsites=[], exceptions=[], miscuts=2):
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
            if cut == 0 and len(frags) == 0:
                frags.append(chop+seq)
            elif cut == 0 and chop and seq:
                frags.append(chop+seq)
            #print(chop)
            elif len(seq) != 0 : frags.append(seq)
            elif len(seq) == 0 and chop : frags.append(chop) #special case for
                                   #if last amino acid is one of the cut sites
            break #no more cuts
        if seq[0] not in exceptions:
            frags.append(chop)
            chop = ''
    merged_list = list()
    no_met = list()
    for k in range(0, miscuts):
        for ix, chunk in enumerate(rolling_window(frags, k+2)):  # +2 to join adjacent
            merged_list.append(''.join(chunk))
            if ix == 0:
                no_met = list()
                for ix, c in enumerate(chunk):
                    if ix == 0:
                        no_met.append(c[1:])
                    else:
                        no_met.append(c)
        merged_list.append(''.join(no_met))
    frags.append(frags[0][1:]) # chop off methionine
    nomiscuts_len = len([  x for x in frags if len(x) >= minlen  ])
    return [ x for x in frags+merged_list if len(x) >= minlen ], nomiscuts_len


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





def seq_modi(sequence, modifications, keeplog=True):

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

    also modifies sequence to fill in Xs with the predicted modification
    example : AAAx(P)R --> AAAPR
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
                        if keeplog:
                            logging.warning('New modification {} that'\
                                            ' is not found in sequence {}'.format(
                                                modi, sequence))

    sequence = ''.join(seqlist)
    if not seqmodi:
        seqmodi = sequence
    return sequence, seqmodi, modi_len, label

def _count_modis_maxquant(modi):
    if modi == 'Unmodified':
        return 0
    return modi.count(',') + 1  # multiple modis are separated by a comma

def count_modis_maxquant(df):
    return df.apply(lambda x: _count_modis_maxquant(x['Modifications']), axis=1)

def calculate_miscuts(seq, targets=None):
    """Calculates number of miscuts for a given sequence
    using amino acids given in targets"""
    # TODO: implement support for cases such as trypsin/P that don't count miscuts with next aa
    # a Proline, for example
    miscuts = sum(seq.count(x) for x in targets)
    if not any(seq[-1] == x for x in 'KR'):  # then at the C terminal
        return miscuts
    return miscuts -1
