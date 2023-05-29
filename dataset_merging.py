import pandas as pd
import numpy as np

def merge_bgc(path = 'data/bgcs/'):
    df_dict = dict()

    df_dict['CDS_motifs'] = pd.read_csv(path + 'all_features_CDS_motifs.csv').iloc[:,2:]
    df_dict['CDS_smCOG'] = pd.read_csv(path + 'all_features_CDS_smCOG.csv').iloc[:,2:]

    df_dict['PFAM_ID'] = pd.read_csv(path + 'all_features_PFAM_ID.csv').iloc[:,2:]

    df_dict['labels'] = pd.read_csv(path + 'prod_activity_identif.csv').iloc[:,1:]

    df_merged = pd.merge(df_dict['CDS_motifs'], df_dict['CDS_smCOG'], on = 'BGC_id', how = 'inner')
    df_merged = pd.merge(df_merged, df_dict['PFAM_ID'], on = 'BGC_id', how = 'inner')
    df_merged = pd.merge(df_merged, df_dict['labels'], on = 'BGC_id', how = 'inner')

    df_merged.insert(0, 'cluster_name', df_merged.pop('cluster_name'))
    df_merged.insert(1, 'BGC_id', df_merged.pop('BGC_id'))

    df_merged.to_csv(path + 'data.csv', index=False)

def merge_non_antifungi(path, label):
    df_dict = {}

    #folder = 'data/acc/'

    df_dict['CDS_motifs'] = pd.read_csv(path + 'all_features_CDS_motifs.csv')
    df_dict['CDS_motifs'] = df_dict['CDS_motifs'].rename(columns = {df_dict['CDS_motifs'].columns[0]: 'index'}).drop(df_dict['CDS_motifs'].iloc[:, 1:6], axis=1)
    df_dict['CDS_smCOG'] = pd.read_csv(path + 'all_features_CDS_smCOG.csv')
    df_dict['CDS_smCOG'] = df_dict['CDS_smCOG'].rename(columns = {df_dict['CDS_smCOG'].columns[0]: 'index'}).drop(df_dict['CDS_smCOG'].columns[[1]], axis=1)
    df_dict['PFAM_desc'] = pd.read_csv(path + 'all_features_PFAM_desc.csv')
    df_dict['PFAM_desc'] = df_dict['PFAM_desc'].rename(columns = {df_dict['PFAM_desc'].columns[0]: 'index'}).drop(df_dict['PFAM_desc'].columns[[1]], axis=1)
    df_dict['PFAM_GO'] = pd.read_csv(path + 'all_features_PFAM_GO.csv')
    df_dict['PFAM_GO'] = df_dict['PFAM_GO'].rename(columns = {df_dict['PFAM_GO'].columns[0]: 'index'}).drop(df_dict['PFAM_GO'].columns[[1]], axis=1)
    df_dict['PFAM_ID'] = pd.read_csv(path + 'all_features_PFAM_ID.csv')
    df_dict['PFAM_ID'] = df_dict['PFAM_ID'].rename(columns = {df_dict['PFAM_ID'].columns[0]: 'index'}).drop(df_dict['PFAM_ID'].columns[[1]], axis=1)

    df_merged = None

    for df in df_dict:
        if df_merged is None:
            df_merged = df_dict[df]
        else:
            df_merged = pd.concat(df_merged, df_dict[df], on = 'index', how = 'inner')

    df_merged = df_merged.drop("index", axis=1)
    df_merged.loc[:,'is_antifung'] = label

    df_merged.to_csv(path + 'data.csv', index=False)
