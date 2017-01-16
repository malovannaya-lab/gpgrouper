"""Original, non vectorized code
Keeping around in case there is an error with new code logic"""
import numpy as np
import pandas as pd

def _get_peptides_for_gene(df,usrdata):
    '''
    Lookup info from userdata. Also sorts based on alphabet.
    '''
    IDquery = df['e2g_GeneID']
    #print(IDquery)  # for debugging
    try :
        matches = usrdata[usrdata.psm_GeneID == IDquery].copy()
        matches = matches.sort_values(by=['psm_PSM_IDG'], ascending=True).drop_duplicates(subset='sequence_lower',
                                                                                keep='first',)
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
    return (set(sorted(matches.sequence_lower)), pepts_str, protcount, uniques,
         protcount_S, uniques_S)

def get_peptides_for_gene(genes_df, temp_df):
    """Get peptide sequence information for each gene"""
    peptides = genes_df.apply(_get_peptides_for_gene, args=(temp_df,), axis=1)
    (genes_df['e2g_PeptideSet'], genes_df['e2g_PeptidePrint'], genes_df['e2g_PeptideCount'],
     genes_df['e2g_PeptideCount_u2g'], genes_df['e2g_PeptideCount_S'],
     genes_df['e2g_PeptideCount_S_u2g']) = list(zip(*peptides))
    return genes_df

def _get_psms_for_gene(gene_df_ID, data, EXPTechRepNo=1):
    total = data[data['psm_GeneID']==gene_df_ID]['psm_PSM_UseFLAG'].sum()/EXPTechRepNo

    total_u2g = (data[(data['psm_GeneID']==gene_df_ID) & (data['psm_GeneCount' ] == 1)]
                 ['psm_PSM_UseFLAG']).sum()/EXPTechRepNo

    total_S = (data[(data['psm_GeneID']==gene_df_ID) & (data['psm_PSM_IDG' ] < 4)]
               ['psm_PSM_UseFLAG']).sum()/EXPTechRepNo

    total_S_u2g = (data[(data['psm_GeneID']==gene_df_ID) & (data['psm_PSM_IDG' ] < 4) &
                        (data['psm_GeneCount' ] == 1)]
                   ['psm_PSM_UseFLAG']).sum()/EXPTechRepNo
    #for match in matches:
    return total, total_u2g, total_S, total_S_u2g

def get_psms_for_gene(genes_df, temp_df,):
    """Get PSMs information for each gene"""
    psm_info = genes_df.apply(lambda x:
                              _get_psms_for_gene(x['e2g_GeneID'],
                                                 temp_df,
                              ),
                              axis=1)
    (genes_df['e2g_PSMs'], genes_df['e2g_PSMs_u2g'],
     genes_df['e2g_PSMs_S'],genes_df['e2g_PSMs_S_u2g']) = list(zip(*psm_info))
    return genes_df


def _assign_IDG(ionscore, qvalue, ion_score_bins):

    low, med, high = ion_score_bins  # default should be 10, 20, 30
    if (ionscore>= high and qvalue <= .01): IDGout = 1
    elif (ionscore >= high and qvalue <= .05): IDGout = 2
    elif (ionscore >= med and qvalue <= .01): IDGout = 3
    elif (ionscore >= med and qvalue <= .05): IDGout = 4
    elif (ionscore >= low and qvalue <= .01): IDGout = 5
    elif (ionscore >= low and qvalue <= .05): IDGout = 6
    elif (ionscore >= 0 and qvalue <= .01): IDGout = 7
    elif (ionscore >= 0 and qvalue <= .05): IDGout = 8
    else: IDGout  = 9
    return IDGout

def assign_IDG(usrdata):
    """Assign IDG bsaed on combination of
    IonScore and q_value"""
    ion_score_bins = usrdata.filtervalues.get('ion_score_bins', (10, 20, 30))
    usrdata.df['psm_PSM_IDG'] = usrdata.df.apply(lambda x:
                                                 _assign_IDG(x['IonScore'],
                                                             x['q_value'],
                                                             ion_score_bins), axis=1)
    return usrdata


# def _peptidome_matcher(seq, metadata, prot):
#     seq = seq.upper()
#     if not isinstance(metadata, tuple):
#         metadata = tuple()
#     if seq in prot:
#         metadata = set(metadata).union(prot[seq])
#     return tuple(metadata)

# def peptidome_matcher(usrdata, ref_dict):
#     """Matches Sequence column with refseq dictionary
#     returns an empty list if there is no match.
#     returns input DataFrame with a metadata column with a list of named tuples"""
#     usrdata['metadatainfo'] = usrdata.apply(lambda x:
#                                             _peptidome_matcher(x['Sequence'],
#                                                                x['metadatainfo'],
#                                                                ref_dict),
#                                             axis=1)
#     return usrdata


def _flag_AUC_PSM(df, d):
    if d['pep'] =='all' : d['pep'] = float('inf')
    if d['idg'] =='all' : d['idg'] = float('inf')
    #if df['IonScore'] == '' and df['q-Value'] == '':
    if df['Charge'] < d['zmin'] or df['Charge'] > d['zmax'] :
         # if charge is empty (nan), this will not be true
        AUC_flag = 0
        PSM_flag = 0
    elif df['psm_SequenceModiCount'] > d['modi'] :
        AUC_flag = 0
        PSM_flag = 0
    elif np.isnan(df['IonScore']) and np.isnan(df['q_value']):
        AUC_flag = 1 #retains WLs Q2up assignments
        PSM_flag = 0
    elif df['IonScore'] < d['ion_score'] :
        AUC_flag = 0
        PSM_flag = 0
    elif df['q_value'] > d['qvalue'] :
        AUC_flag = 0
        PSM_flag = 0
    elif df['PEP'] > d['pep'] :
        AUC_flag = 0
        PSM_flag = 0
    elif df['psm_PSM_IDG'] > d['idg'] :
        AUC_flag = 0
        PSM_flag = 0
    elif df['psm_Peak_UseFLAG'] == 0 :
        AUC_flag = 0
        #if 'PSMAmbiguity' in df.columns:  # will not work like this
        if df['PSMAmbiguity'].lower() == 'unambiguous':
            PSM_flag = 1
        else:
            PSM_flag = 0
    else:
        AUC_flag = 1
        PSM_flag = 1

    if df['AUC_reflagger'] == 0:
        AUC_flag = 0

    return AUC_flag, PSM_flag

def flag_AUC_PSM(usrdata, filtervalues):
    """Apply AUC and PSM flags per row"""
    flags = usrdata.apply(_flag_AUC_PSM, args=(filtervalues,), axis=1)
    usrdata['psm_AUC_UseFLAG'], usrdata['psm_PSM_UseFLAG'] = list(zip(*flags))
    return usrdata

# def _gene_taxon_map(gene, d):
#     return d.get(gene)

# def gene_taxon_map(usrdata, gene_taxon_dict):
#     """make 'gene_taxon_map' column per row which displays taxon for given gene"""
#     usrdata['psm_TaxonID'] = usrdata.apply(lambda x : _gene_taxon_map(
#         ['psm_GeneID'], gene_taxon_dict), axis=1)
#     return usrdata


def _get_gene_capacity(geneid, database):
    return database[database.GeneID == geneid].capacity.mean()

def get_gene_capacity(genes_df, database, col='e2g_GeneID'):
    """Get gene capcaity from the stored metadata"""
    genes_df['e2g_GeneCapacity'] = genes_df.apply(lambda x:
                                                  _get_gene_capacity(
                                                      x[col],
                                                      database),
                                                  axis=1)
    return genes_df


def _calculate_protein_area(gene_df, usrdata, area_col, normalization, EXPTechRepNo=1):
    """Area column is psm_SequenceArea"""

    matches  = usrdata[(usrdata['psm_GeneID'] == gene_df['e2g_GeneID']) &
                           (usrdata['psm_AUC_UseFLAG']==1)] [
                                [area_col, 'psm_GeneCount',
                                 'psm_PSM_IDG', 'MissedCleavages']]

    uniq_matches = matches[matches['psm_GeneCount']==1]
    uniq_matches_0 = uniq_matches[uniq_matches['MissedCleavages']==0]
    matches_strict = matches[matches['psm_PSM_IDG'] < 4]
    values_max = matches[area_col].sum()
    values_adj = (matches[area_col]/matches['psm_GeneCount']).sum()
    uniq_values_adj = (uniq_matches[area_col]/matches['psm_GeneCount']).sum()
    uniq_values_adj_0 = (uniq_matches_0[area_col]/matches['psm_GeneCount']).sum()
    result = (values_max/normalization, values_adj/normalization,
              uniq_values_adj_0/normalization, uniq_values_adj/normalization)
    return result


def calculate_protein_area(genes_df, temp_df, area_col, normalize):
    """Calculate the area of each protein"""
    areas = genes_df.apply(_calculate_protein_area,
                                    args=(temp_df,
                                           area_col, normalize),
                                     axis=1)
    (genes_df['e2g_nGPArea_Sum_max'],genes_df['e2g_nGPArea_Sum_cgpAdj'],
     genes_df['e2g_nGPArea_Sum_u2g'], genes_df['e2g_nGPArea_Sum_u2g_all']) = list(zip(*areas))
    return genes_df

def _calculate_gene_dstrarea(genes_df, temp_df, normalization):
    gene = genes_df['e2g_GeneID']
    if genes_df['e2g_IDSet'] in [1,2]:
        return (temp_df[(temp_df['psm_GeneID'] == gene) &
                       (temp_df['psm_AUC_UseFLAG'] == 1)]
                .psm_PrecursorArea_dstrAdj).sum()/normalization
    elif genes_df['e2g_IDSet'] == 3:
        return 0

def calculate_gene_dstrarea(genes_df, temp_df, normalize):
    """Calculate distributed area for each gene product"""
    genes_df['e2g_nGPArea_Sum_dstrAdj'] = genes_df.apply(_calculate_gene_dstrarea,
                                                      args=(temp_df,
                                                            normalize,),
                                                      axis=1)
    return genes_df
