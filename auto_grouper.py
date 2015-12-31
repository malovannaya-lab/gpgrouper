import os
import time
import threading
import argparse
#from collections import defaultdict
from configparser import ConfigParser
#from datetime import datetime
import pandas as pd
#import database_config as db
import pygrouper
import bcmproteomics as bcm

def experiment_checker():
    """Looks for experiments that have records in ispec but have not been grouped yet
    """
    conn = bcm.filedb_connect()
    print('Looking for new experiments from iSPEC...')
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

    sql = ("SELECT {} FROM iSPEC_BCM.iSPEC_ExperimentRuns "
           "WHERE exprun_Grouper_EndFLAG != 1 AND "
           "exprun_Grouper_StartFLAG != 1 AND "
           "exprun_Grouper_FailedFLAG != 1").format(', '.join(exprun_cols))
    #sql = 'SELECT {} FROM iSPEC_BCM.iSPEC_ExperimentRuns '\
    #      'WHERE exprun_EXPRecNo in ({})'.format(', '.join(exprun_cols),
    #                                             ', '.join(result))
    rundata = pd.read_sql(sql, conn, index_col='exprun_EXPRecNo')
    if len(rundata) == 0:
        return rundata
    rundata = rundata[~rundata['exprun_EXPRunNo'].isnull()]  # only keep if has run number
    sql = ("SELECT exp_EXPRecNO, exp_Digest_Type, exp_Digest_Enzyme, exp_CreationTS, "
           "exp_AddedBy from iSPEC_BCM.iSPEC_Experiments "
           "WHERE exp_EXPRecNo in ({})").format(', '.join([str(rec) for
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

    return rundata
    #sql = ('SELECT record_no, run_no, search_no FROM experimentruns WHERE '
    #       'record_no in ({})').format(', '.join([str(rec) for rec in rundata.index.tolist()]))
    #local_conn = db.get_connection()
    #local_recs = pd.read_sql(sql, local_conn, index_col='record_no')
    #newexps = defaultdict(list)
    #for ix, row in rundata.iterrows():
    #    ix = int(ix)
    #    if isinstance(row.exp_CreationTS, str):
    #        try:
    #            creation = datetime.strptime(row.exp_CreationTS, '%Y/%m/%d %H:%M:%S')
    #        except ValueError:
    #            creation = None
    #    else:
    #        creation = row.exp_CreationTS
    #    in_local_db = False
    #    if len(local_recs[(local_recs.index==ix) & (local_recs.run_no == row['exprun_EXPRunNo']) &
    #                      (local_recs.search_no==row['exprun_EXPSearchNo'])]) == 0:
    #        newexps[ix].append(
    #            {'run_no':row.exprun_EXPRunNo,
    #             'search_no': row.exprun_EXPSearchNo, 'taxon': row.exprun_TaxonID,
    #             'addedby': row.exprun_AddedBy, 'creation_ts': creation,
    #             'purpose': row.exprun_Purpose, 'label': row.exprun_LabelType,
    #             'techrep': row.exprun_nTechRepeats, 'instrument':row.exprun_MS_Instrument,
    #             'msexperimenter': row.exprun_MS_Experimenter, 'mscomment':row.exprun_MS_Experimenter,
    #             'quant': row.exprun_Search_QuantSource, 'searchexperimenter': row.exprun_Search_Experimenter}
    #            )

    #if newexps:
    #    print('Updating experiment records')
    #    db.add_experiments(newexps)
    #conn.close()

def file_checker(INPUT_DIR, OUTPUT_DIR, maxqueue):
    """Docstring
    """
    validfiles = [f for f in os.listdir(INPUT_DIR) if '.txt' in f]
    setups, usrfiles = [], []
    #ungrouped = db.get_ungrouped_exps()
    #session = db.make_session()
    ungrouped = experiment_checker()  # may be empty, is ok
    usrfilesize = 0  # keep track of usrfile size so we don't group too many at once
    MAX_SIZE = 2796774193548.3867
    queue_size = 0
    for recno, exp in ungrouped.iterrows():  # the index is the record number
        recno = int(recno)  # int to remove decimal to match file name
        #try:
        #    query = session.query(db.Experiment).filter_by(record_no=exp.record_no).one()
        #    added_by = query.added_by
        #except NoResultFound:  # should not occur do to how database is set up
        #    added_by = ''

        setup = {'EXPRecNo': recno,
                 'EXPRunNo': exp.exprun_EXPRunNo,
                 'EXPSearchNo': exp.exprun_EXPSearchNo,
                 'taxonID': exp.exprun_TaxonID,
                 'EXPQuantSource': exp.exprun_Search_QuantSource,
                 'AddedBy': exp.exprun_AddedBy,
                 'EXPTechRepNo': exp.exprun_nTechRepeats,
                 'EXPLabelType': exp.exprun_LabelType,
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
        if setup and usrfile and (usrfilesize <= MAX_SIZE) and queue_size<=maxqueue:  # if we have both,
            setups.append(setup)                     # a cap on max files to group at once
            usrfiles.append(usrfile)
            queue_size += 1
            print('Found experiment {} from datafile'\
                  ' {}.'.format(setup['EXPRecNo'], usrfile))
            conn = bcm.filedb_connect()
            cursor = conn.cursor()
            cursor.execute("""UPDATE iSPEC_ExperimentRuns 
            SET exprun_Grouper_StartFLAG = ?
            WHERE exprun_EXPRecNo= ? AND 
            exprun_EXPRunNo = ? AND 
            exprun_EXPSearchNo = ?
            """, 1, setup['EXPRecNo'],
                           setup['EXPRunNo'], setup['EXPSearchNo'])
            conn.commit()
    if len(usrfiles) > 0:
        pygrouper.main(usrfiles=usrfiles, exp_setups=setups, automated=True,
                       inputdir=INPUT_DIR, outputdir=OUTPUT_DIR, usedb=True)
    #session.close()

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
    print('parser time first')
    parser = argparse.ArgumentParser(
        description="""Script to run pygrouper automatically.
        Requires access to iSPEC (via pyodbc) and an iSPEC login.""",

        epilog=""" """
        )

    parser.add_argument('-m','--max', type=int, default=999,
                        help=('Set a maximum number of experiments to queue '
                              'at once. Useful if you want to run pygrouper '
                              'in small batches on multiple computers.')
                        )

    args = parser.parse_args()

    configparser = ConfigParser()
    configparser.read('py_config.ini')
    INPUT_DIR = configparser['directories']['inputdir']  # where to look
                                                   # for files to group
    OUTPUT_DIR = configparser['directories']['outputdir']

    while True:
        INTERVAL = 60 * 60
        schedule(INTERVAL, [INPUT_DIR, OUTPUT_DIR, args.max])
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
