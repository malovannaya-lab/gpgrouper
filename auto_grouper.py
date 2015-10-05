import os
import time
import sys
import threading
from configparser import SafeConfigParser

import database_config as db
import pygrouper

ispecf = '4PyGrouper_ExpRunDump.xlsx'

def file_checker(INPUT_DIR, OUTPUT_DIR):
    validfiles = [f for f in os.listdir(INPUT_DIR) if '.txt' in f]
    setups, usrfiles = [], []
    ungrouped = db.get_ungrouped_exps()
    session = db.make_session()
    usrfilesize = 0  # keep track of usrfile size so we don't group too many at once
    MAX_SIZE = 2796774193548.3867
    for exp in ungrouped:
        try:
            query = session.query(db.Experiment).filter_by(record_no=exp.record_no).one()
            added_by = query.added_by
        except NoResultFound:  # should not occur do to how database is set up
            added_by = ''

        setup = {'EXPRecNo': exp.record_no,
                     'EXPRunNo': exp.run_no,
                     'EXPSearchNo': exp.search_no,
                     'taxonID': exp.taxonid,
                     'EXPQuantSource': exp.quant_source,
                     'AddedBy': added_by,
                     'EXPTechRepNo': exp.tech_repeat,
                     'EXPLabelType': exp.label_type
                     }
        expfilematch = str(setup['EXPRecNo'])+'_'+str(setup['EXPRunNo'])
        usrfilelist = [f for f in validfiles if f.startswith(expfilematch)]
        if len(usrfilelist) == 1: # ensure we have just one match
            usrfile = usrfilelist[0]
        elif len(usrfilelist) > 1:
            usrfilelist = [f for f in usrfilelist if 'all' in f]  # 'all' should have all fractions combined

            if len(usrfilelist) == 1:
                usrfile == usrfilelist[0]
            else:
                usrfile = None
                print('More than one match for experiment {}'\
                      ', skipping'.format(expfilematch))
        elif len(usrfilelist) == 0: # file not here yet
            usrfile = None
        if usrfile:
            usrfilesize += os.stat(os.path.join(INPUT_DIR, usrfile)).st_size
            

        if setup and usrfile and (usrfilesize <= MAX_SIZE):  # if we have both, 
            setups.append(setup)                     # a cap on max files to group at once
            usrfiles.append(usrfile)
            print('Found experiment {} from datafile'\
                  ' {}.'.format(setup['EXPRecNo'], usrfile))
    if len(usrfiles) > 0:
        pygrouper.main(usrfiles=usrfiles, exp_setups=setups, automated=True,
                       inputdir=INPUT_DIR, outputdir=OUTPUT_DIR, usedb=True)
    session.close()
        
def schedule(INTERVAL, args):
    print('{} : Checking for new experiments.'.format(time.ctime()))
    INPUT_DIR = args[0]
    db.get_ispec_experiment_info(os.path.join(INPUT_DIR,ispecf), todb=True, norepeats=True)
    file_checker(*args)
    global thread
    thread = threading.Timer(INTERVAL, schedule, [INTERVAL, args])
    print('{} : Sleeping... Press Enter to wake or [exit] to'\
          ' end.'.format(time.ctime()), end='\n\n')
    thread.start()
    

if __name__ == '__main__':
    parser = SafeConfigParser()
    parser.read('py_config.ini')
    INPUT_DIR = parser['directories']['inputdir']  # where to look
                                                   # for files to group
    OUTPUT_DIR = parser['directories']['outputdir']
    

    while True:
        INTERVAL = 60 * 60
        schedule(INTERVAL, [INPUT_DIR, OUTPUT_DIR])
        #print('{} : Sleeping... Press Enter to wake or [exit] to
        #end.'.format(time.ctime()))
        usr = input()  # break from schedule interval and manually call
        if usr.lower() == 'exit':
            print('Goodbye.\n')
            thread.cancel()
            break
        else:
            thread.cancel()  # Stop the timer and cancel the execution of the
            #timer's action


    
