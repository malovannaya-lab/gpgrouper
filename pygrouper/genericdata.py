#!/usr/bin/env python
#Purpose: Download and reformat generic protein- and gene-related databases from NCBI and EBI.
#Laboratory: Malovannaya Lab, Baylor College of Medicine, One Baylor Plaza BCM 125, Houston, TX 77025, USA
#Author(s): Bhoomi Bhatt; Alexander Saltzman (edit); Anna Malovannaya (edit)
#Date: 20150806

import os
import sys
import shutil
import ftplib
from ftplib import FTP
import gzip
import datetime
import re
import csv


now = datetime.datetime.now()
dt = now.strftime('%Y%m%d')

org_dict = {9606: 'HUMAN',
            10090: 'MOUSE'}
ftp_dict = {9606: 'refseq/H_sapiens/mRNA_Prot',
            10090: 'refseq/M_musculus/mRNA_Prot'}
def cleanup(path):
    """Cleanup temp directory and path"""
    shutil.rmtree(path)

def download_ebi_files(path='.', taxa=None):
    """Download necessary files """
    #EBI downloads: (1) IPREntry (entry.list); (2) protein2ipr; (3) idmapping.

    #Downloading IPRentry from EBI
    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login()
    ftp.cwd('/pub/databases/interpro')
    out = os.path.join(path, 'entry_.tab')
    ftp.retrbinary('RETR entry.list', open(out, 'wb').write)

    #Downloading protein2ipr from EBI
    out = os.path.join(path, 'protein2ipr_.gz')
    if not os.path.exists(out):
        ftp.retrbinary('RETR protein2ipr.dat.gz', open(out, 'wb').write)

    #Downloading idmapping from EBI
    ftp.cwd('/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism')
    for taxon in taxa:
        out = os.path.join(path,
                           '{}_{}_idmapping_selected.gz'.format(org_dict.get(taxon), taxon))

        if not os.path.exists(out):
            ftp.retrbinary('RETR {}_{}_idmapping_selected.tab.gz'.format(org_dict.get(taxon), taxon),
                           open(out, 'wb').write)

    ftp.quit()

def download_ncbi_files(path='.', taxa=None):
    """Download necessary files"""
    #NCBI downloads: (1) RefSeq (*protein.faa.gz); (2) homologene.data; (3) gene directory files.
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    filematch = '*faa.gz'
    def get_taxa_data():
        for filename in ftp.nlst(filematch):
            fname, extension = os.path.splitext(filename)
            out = os.path.join(path,
                               '{}{}'.format(fname, extension))
            if not os.path.exists(out):
                fhandle = open(out, 'wb')
                print ('Getting files ' + filename)
                ftp.retrbinary('RETR {}'.format(filename),  fhandle.write)
                fhandle.close()

    for taxon in taxa:
        ftp.cwd(ftp_dict.get(taxon))
        get_taxa_data()

    #Downloading homologene from NCBI homologene FTP
    ftp.cwd('/pub/HomoloGene')  #get latest built
    latest = sorted([x for x in ftp.nlst() if x.startswith('build')], reverse=True)[0]
    ftp.cwd('/pub/HomoloGene/{}'.format(latest))  #get latest built
    out = os.path.join(path,
                       'homologene_download.tab')

    if not os.path.exists(out):
        ftp.retrbinary('RETR homologene.data',
                    open(out, 'wb').write)
    ftp.cwd('/gene/DATA')
    filematch = '*.gz'
    for filename in ftp.nlst(filematch):
        fname, extension= os.path.splitext(filename)
        out = os.path.join(path,
                           '{}{}'.format(fname, extension))
        if not os.path.exists(out):
            fhandle = open(out, 'wb')
            print('Getting ' + filename)
            ftp.retrbinary('RETR {}'.format(filename),  fhandle.write)
            fhandle.close()
    for outpath, retbin in [(os.path.join(path, 'go_process.dtd'), 'RETR go_process.dtd'),
                         (os.path.join(path, 'go_process.xml'), 'RETR go_process.xml'),
                         (os.path.join(path, 'mim2gene_medgen'), 'RETR mim2gene_medgen')]:
        if not os.path.exists(outpath):
            ftp.retrbinary(retbin, open(outpath).write)


# Unzip Files:
def unzip(inputf, outputf):
    "unzip a file"
    with gz.open(inputf, 'rb') as zipf, open(outputf, 'wb') as outf:
        for line in zipf:
            outf.write(line)

def unzip_all(path='.'):
    """unzip all gz files in a given directory"""
    for entry in os.scandir(path):
        if entry.is_file() and entry.name.endswith('.gz'):
            entry_base, _ = os.path.splitext(entry.name)
            unzip(entry.name, entry_base + '.tab')

# Reformat Files:
def entrylist_formatter(filein=None, fileout=None, path='.'):

    if filein is None:
        filein = 'entry.tab'
    if fileout is None:
        fileout = 'IPRentry.tab'

    category = ''
    with open(filein, 'r') as fin, open(fileout, 'w') as fout:
        for line in fin:
            if line.startswith('IPR'):
                linesplit = [line[0:10], line[10:]]
                fout.write('{}\t{}'.format(category,
                                           '\t'.join(linesplit[0],
                                                     linesplit[1])))
            else: # start of a new category
                category = line.strip()

def protein2ipr_formatter(filein=None, fileout=None, path='.'):

    if filein is None:
        filein = 'protein2ipr.tab'
    if fileout is None:
        fileout = 'protein2ipr_uni2ipr.tab'
    filein = os.path.join(path, filein)
    fileout = os.path.join(path, fileout)
    with open(filein, 'r') as fin, open(fileout, 'w') as fout:
        for line in fin:
            linesplit = line.split('\t')
            fout.write('\t'.join(linesplit[0:2]))  # only IPRIDs and description

def idmapping_formatter(taxonid, fileout=None, path='.'):

    if fileout is None:
        fileout = 'idmapping_{}.tab'.format(taxonid)
    fileout = os.path.join(path, fileout)

    inputfiles = list()
    searchterm = org_dict.get(taxonid)
    filepat = re.compile('{}_{}_idmapping_selected.tab'.format(searchterm,
                                                                  taxonid,),
                         re.IGNORECASE)

    heading=['UniProt_AC', 'UniProtKB_ID', 'GeneID', 'RefSeq', 'GI',
             'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR',
             'taxonID', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL-CDS', 'Ensembl',
             'Ensembl_TRS', 'Ensembl_PRO', 'Additional PubMed']
    # Output file header line

    for entry in os.scandir(path):
        if entry.is_file and filepat.search(entry):
            inputfiles.append(entry.name)
    with open(fileout, 'w') as fout:
        fout.write('\t'.join(heading))
        fout.write('\n')
        for inputf in inputfiles:
            with open(inputf, 'r') as fin:
                for line in fin:
                    fout.write(line)

def gene_info_formatter(taxonid, filein=None, fileout=None):

    if filein is None:
        filein = 'gene_info.tab'
    if fileout is None:
        fileout = 'gene_info_formatted_{}.tab'.format(taxonid)
    filein = os.path.join(path, filein)
    fileout = os.path.join(path, fileout)

    heading=['TaxonID', 'GeneID', 'GeneSymbol', 'LocusTag', 'GeneSynonyms', 'dbXrefs',
             'GeneChrLocation', 'GeneMapLocation', 'GeneDescription', 'GeneType',
             'GeneSymbolFromNomenclatureAuthority', 'GeneFullNameFromNomenclatureAuthority',
             'GeneNomenclatureStatus', 'GeneOtherDesignations', 'GeneModificationDate']

    with open(filein, 'r') as fin, open(fileout, 'w') as fout:
        fout.write('\t'.join(heading))
        fout.write('\n')
        for line in fin:
            linesplit = line.strip().split('\t')
            if linesplit[0] == taxonid:
                fout.write('\t'.join(linesplit[0:16]))
                fout.write('\n')

def append_all_files(taxon, path='.'):
    """get a list of files that match a given searchterm"""
    searchterm = org_dict.get(taxon)
    inputfiles = list()
    searchterm = re.compile('{}\.[0-9=+\.protein.faa.tab'.format(searchterm),
                            re.IGNORECASE)
    for entry in os.scandir(path):
        if entry.is_file and searchterm.search(entry.name):
            inputfiles.append(os.path.append(path, entry.name))
    return inputfiles


def refseq_dict(inputfiles):
    """build a dict with refseq description => a lit with
    proteingi, refseq accession, and then fasta description"""

    refseq = {}
    for ifilename in inputfiles:
        fasta_parser = SeqIO.parse(ifilename, "fasta")
        with open(ifilename,'r') as f:
            for sublst in fasta_parser:
                description_split=sublst.description.split('|')
                refseq[description_split[4].strip()]='\t'.join([description_split[1],
                                                                description_split[3],
                                                                description_split[4],
                                                                str(sublst.seq)])
    return refseq


def gene2accession_formatter(taxonid, filein=None, fileout=None, path='.'):

    if filein is None:
        filein = 'gene2accession.tab'
    if fileout is None:
        fileout = 'gene2accession_formatted_{}.tab'.format(taxonid)
    filein = os.path.join(path, filein)
    fileout = os.path.join(path, fileout)
    with open(filein,'r') as filein, open(fileout,'w') as fileout:
        fileout.write("TaxonID\tGeneID\tProteinGI\n") #Add tabs between the headers
        for line in filein:
            linesplit = line.strip().split('\t')
            if linesplit[0] == taxonid:
                fileout.write('\t'.join([linesplit[0],linesplit[1],linesplit[6]]))
                fileout.write('\n')

def g2acc_dict(taxonid, filein=None, path='.'):

    if filein is None:
        filein = 'gene2accession_formatted_{}.tab'.format(taxonid)
    filein = os.path.join(filein)
    g2a = dict()
    with open(filein, 'r') as f:
        for line in f:
            linesplit=line.strip().split('\t')
            g2a[linesplit[2].strip()]=linesplit[1]   #keys=proteingi and values=geneid
    return(g2a)

def homologene_formatter(taxonid, filein=None, fileout=None, path='.'):
    'just adds tabs to a header'
    if filein is None:
        filein = 'homologene_download.tab'
    filein = os.path.join(path, filein)
    header = ['hid_HomologeneID', 'hid_TaxonID', 'hid_GeneID', 'hid_GeneSymbol', 'hid_ProteinGI',
              'hid_ProteinAccession']
    with open(filein, 'r') as fin, open(fileout, 'w') as fout:
            fout.write('\t'.join(header))
            fout.write('\n')
            for line in fin:
                fout.write(line)

def hid_dict(filein=None, path='.'):
    if filein is None:
        filein = 'homologene_download.tab'
    filein = os.path.join(path, filein)
    hid = dict()
    with open (filein, 'r') as f:
        for line in f:
            linesplit=line.strip().split('\t')
            hid[linesplit[4].strip()] = linesplit[0] #keys=proteingi and values=homologeneid
    return hid


def file_input(inputfiles=None, refseq=None, g2a=None, hid=None, taxonid=None):
    """writing the reseq, homologene ids and gene id to a set variable to remove redundant lines"""
    if any(x is None for x in (inputfiles, refseq, g2a, hid, taxonid)):
        raise TypeError("one or few input missing")
        return

    lines_seen=set()
    for ifilename in inputfiles:
        with open(ifilename,'r') as f:
            for line in f:
                if line.startswith('>'):
                    linesplit = line.strip().split('|')
                    try:
                        if linesplit[1].strip() in g2a and linesplit[1].strip() in hid and linesplit[4].strip() in refseq: #Same lines are present in different files and will match to the dictionary; making a set will prevent doing so
                            lines_seen.add('\t'.join([TaxonID,g2a[linesplit[1].strip()],hid[linesplit[1].strip()], refseq[linesplit[4].strip()]]))
                        else:
                            lines_seen.add('\t'.join([TaxonID,g2a[linesplit[1].strip()],'', refseq[linesplit[4].strip()]]))
                    except KeyError:
                        pass
        f.close()
    return lines_seen

def file_write(taxon, lines_seen, fileout=None, path='.'):
    if fileout is None:
        fileout = '{}_{}_refseq_forgrouper.tab'.format(taxon, dt)
    fileout = os.path.join(path, fileout)
    heading=['TaxonID','GeneID','HomologeneID','ProteinGI','Refseq_Acc','Description','Fasta']
    with open(fileout, 'r'):
        fileout.write('\t'.join(heading))
        fileout.write('\n')
        for line in lines_seen:
            fileout.write(line)
            fileout.write('\n')

def refseq_formatter_4_mascot(taxon, filein=None, fileout=None):
    """write the refseq file to another file compatible with mascot search"""

    if filein is None:
        filein = '{}_{}_refseq_forgrouper.tab'.format(taxon, dt)
    if fileout is None:
        fileout = '{}_{}_refseq_forMascot.tab'.format(taxon, dt)
    filein = os.path.join(path, filein)
    fileout = os.path.join(path, fileout)
    with open(filein, 'r') as fin, open(fileout, 'w') as fout:
        next(f)
        for line in f:
            linesplt = line.split('\t')
            fout.write(">gi|{}|ref|{}|{}\n".format('linesplt[3]',linesplt[4],linesplt[5]))
            fout.write(linesplt[6])
