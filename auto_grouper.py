import os
import time
import threading
from collections import defaultdict
from configparser import ConfigParser
from datetime import datetime
import pandas as pd
import database_config as db
import pygrouper
import bcmproteomics as bcm

#ispecf = '4PyGrouper_ExpRunDump.xlsx'

def experiment_checker():
    """Looks for experiments that have records in ispec but have not been grouped yet
    """
    conn = bcm.filedb_connect()
    print('Looking for ne experiments from iSPEC...')
    #sql = 'SELECT exprun_EXPRecNo from iSPEC_BCM.iSPEC_ExperimentRuns where exprun_cGeneCount=0'
    #cursor = conn.execute(sql)
    #result = cursor.fetchall()
    #if result:
    #    result = [str(r) for recno in result for r in recno] # unpack the nested list
    #elif len(result)==0:
    #    return
    exprun_cols = ['exprun_EXPRecNo', 'exprun_EXPRunNo', 'exprun_EXPSearchNo',
                   'exprun_TaxonID', 'exprun_AddedBy', 'exprun_LabelType',
                   'exprun_nTechRepeats',
                   'exprun_Search_QuantSource', 'exprun_Purpose',
                   'exprun_MS_Instrument', 'exprun_MS_Experimenter',
                   'exprun_Search_Experimenter', ]
    sql = 'SELECT {} FROM iSPEC_BCM.iSPEC_ExperimentRuns '\
          'WHERE exprun_cGeneCount=0'.format(', '.join(exprun_cols),)
    #sql = 'SELECT {} FROM iSPEC_BCM.iSPEC_ExperimentRuns '\
    #      'WHERE exprun_EXPRecNo in ({})'.format(', '.join(exprun_cols),
    #                                             ', '.join(result))
    rundata = pd.read_sql(sql, conn, index_col='exprun_EXPRecNo')
    rundata = rundata[~rundata['exprun_EXPRunNo'].isnull()]  # only keep if has run number
    sql = ('SELECT exp_EXPRecNO, exp_Digest_Type, exp_Digest_Enzyme, exp_CreationTS '
           'from iSPEC_BCM.iSPEC_Experiments '
           'WHERE exp_EXPRecNo in ({})').format(', '.join([str(rec) for
                                                           rec in rundata.index.tolist()]))
    recdata = pd.read_sql(sql, conn, index_col='exp_EXPRecNo')
    rundata = rundata.join(recdata)  # join on index, the record number

    intfields = ['exprun_EXPRunNo', 'exprun_EXPSearchNo',
                 'exprun_TaxonID', 'exprun_nTechRepeats']
    for field in intfields:
        rundata[field] = rundata[field].fillna(0)
        rundata[field] = rundata[field].astype('int')  # get rid of decimals

    for col in [c for c in rundata.columns if c not in intfields]:
        rundata[col] = rundata[col].fillna('')
    sql = ('SELECT record_no, run_no, search_no FROM experimentruns WHERE '
           'record_no in ({})').format(', '.join([str(rec) for rec in rundata.index.tolist()]))
    local_conn = db.get_connection()
    local_recs = pd.read_sql(sql, local_conn, index_col='record_no')
    newexps = defaultdict(list)
    for ix, row in rundata.iterrows():
        ix = int(ix)
        if isinstance(row.exp_CreationTS, str):
            try:
                creation = datetime.strptime(row.exp_CreationTS, '%m/%d/%Y %H:%M:%S')
            except ValueError:
                creation = None
        else:
             creation = row.exp_CreationTS
        in_local_db = False
        if len(local_recs[(local_recs.index==ix) & (local_recs.run_no==row['exprun_EXPRunNo']) &
                          (local_recs.search_no==row['exprun_EXPSearchNo'])]) == 0:
            newexps[ix].append(
                {'run_no':row.exprun_EXPRunNo,
                 'search_no': row.exprun_EXPSearchNo, 'taxon': row.exprun_TaxonID,
                 'addedby': row.exprun_AddedBy, 'creation_ts': creation,
                 'purpose': row.exprun_Purpose, 'label': row.exprun_LabelType,
                 'techrep': row.exprun_nTechRepeats, 'instrument':row.exprun_MS_Instrument,
                 'msexperimenter': row.exprun_MS_Experimenter, 'mscomment':row.exprun_MS_Experimenter,
                 'quant': row.exprun_Search_QuantSource, 'searchexperimenter': row.exprun_Search_Experimenter}
                )

    if newexps:
        print('Updating experiment records')
        db.add_experiments(newexps)
    conn.close()

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
                 'EXPLabelType': exp.label_type,
                 }
        expfilematch = str(setup['EXPRecNo'])+'_'+str(setup['EXPRunNo'])+'_'
        usrfilelist = [f for f in validfiles if f.startswith(expfilematch)]
        if len(usrfilelist) == 1: # ensure we have just one match
            usrfile = usrfilelist[0]
        elif len(usrfilelist) > 1:
            usrfilelist = [f for f in usrfilelist if 'all' in f]
            # 'all' should have all fractions combined

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
    experiment_checker()
    #db.get_ispec_experiment_info(os.path.join(INPUT_DIR,ispecf), todb=True, norepeats=True)
    file_checker(*args)
    global thread
    thread = threading.Timer(INTERVAL, schedule, [INTERVAL, args])
    print('{} : Sleeping... Press Enter to wake or [exit] to'\
          ' end.'.format(time.ctime()), end='\n\n')
    thread.start()

if __name__ == '__main__':
    parser = ConfigParser()
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
