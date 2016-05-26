import sys
import os
import shutil
import argparse
from itertools import repeat
# from pygrouper import grouper
from .. import pygrouper
from pygrouper.containers import UserData
#sys.path.append('..')  # access to pygrouper
#import pygrouper


#print('Initializing pygrouper testing...')
#if not os.path.isdir('./testresults'):
#    print('Removing previous test results.')
#    shutil.rmtree('./testresults')
#    os.mkdir('./testresults')

BASE_DIR = os.path.dirname(os.path.realpath(__file__))

def grab_one():
    setups = [{'EXPRecNo': 30595,
               'EXPRunNo': 1,
               'EXPSearchNo': 2,
               'taxonID': 9606,
               'EXPQuantSource': 'AUC',
               'AddedBy': 'test',
               'EXPTechRepNo': 1,
               'EXPLabelType': 'none'
    },
    ]
    files = ['30595_1_2_3T3_LM2_5_5_PROF_75MIN_all.txt']
    files = [os.path.join(BASE_DIR, f) for f in files]
    return (files, setups)

def get_quick_data():
    data = UserData('30490_1_EQP_6KiP_all.txt')
    data.recno = 30490
    data.taxonid = 9606
    data.quant_source = 'AUC'
    data.added_by = 'test'
    data.labeltype = 'none'
    data.usedb = False
    data.indir = BASE_DIR
    data.outdir = './testresults'

    return data

def get_prof_data():
    data = UserData('30404_1_QEP_ML262_75min_020_RPall.txt')
    data.recno = 30404
    data.taxonid = 20150811
    data.quant_source = 'AUC'
    data.added_by = 'test'
    data.labeltype = 'none'
    data.indir = BASE_DIR
    data.outdir = './testresults'

    return data

def grab_data(quick, prof, tmt):
    quick_data = get_quick_data()
    tmt_data = get_quick_data()
    tmt_data.labeltype = 'TMT'
    tmt_data.datafile = '30490_1_EQP_6KiP_all_TMT.txt'
    prof_data = get_prof_data()

    if quick:
        return [quick_data]
        # files = ['30490_1_EQP_6KiP_all.txt']
    elif tmt:
        return [tmt_data]
    elif prof:
        return [prof_data]
    else:
        return [quick_data, tmt_data, prof_data]

def runtest(quick=False, prof=False, tmt=False, inputdir=None, **kwargs):
    inputdir = BASE_DIR
    usrdatas = grab_data(quick, prof, tmt)
    if inputdir is None:
        inputdir = os.getcwd()
    print('Files for testing :')
    files = [usrdata.datafile for usrdata in usrdatas]
    for f in files:
        print(f)
    print('running test')
    pygrouper.main(usrdatas,
                   usedb=False, **kwargs)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""Script to test pygrouper""",

        epilog=""" """
        )

    parser.add_argument('-q','--quick', action='store_true',
                        help='Run one small file.')
    parser.add_argument('-p','--profile', action='store_true',
                        help='Run one large file.')
    parser.add_argument('-t', '--test', action='store_true',
                        help='Run 1 specific test file.')

    args = parser.parse_args()
    runtest(quick=args.quick, prof=args.profile, testone=args.test)
