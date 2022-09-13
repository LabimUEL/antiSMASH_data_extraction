import pandas as pd
import numpy as np

df_dict = dict()

df_dict['CDS_motifs'] = pd.read_csv('data/all_features_CDS_motifs.csv').iloc[:,2:]
df_dict['CDS_smCOG'] = pd.read_csv('data/all_features_CDS_smCOG.csv').iloc[:,2:]

pfam_df = pd.read_csv('data/all_features_PFAM_ID.csv').iloc[:,2:]

all_GO = set()

with open('data/pfam2go.txt', 'r') as pfam2go:
    
    pfam2g0_dict = {}
    
    for line in pfam2go:
        if line.startswith('!') == False:
            pfam = line.split(' > ')[0].split(' ')[0].split(':')[1]
            GO_term = line.split(' > ')[1].split(' ; ')[1].strip()   
            all_GO.add(GO_term)         
            if pfam not in pfam2g0_dict:
                pfam2g0_dict[pfam] = {GO_term}
            else:
                pfam2g0_dict[pfam].add(GO_term)

pfam_cols = list(pfam_df.columns[['PF' in col for col in pfam_df.columns]])

go_df = pd.DataFrame(np.zeros((len(pfam_df), len(all_GO))), columns = all_GO)

for i in range(len(pfam_df)):
    for pfam in pfam_cols:
        if pfam_df.at[i, pfam] > 0:

            pfam_str = pfam.split('.')[0]

            if pfam_str in pfam2g0_dict:

                for GO_term in pfam2g0_dict[pfam_str]:
                    go_df.at[i, GO_term] += 1

go_df = go_df#.loc[:, (go_df != 0).any(axis=0)] # Remove GO terms not found in any cluster

df_dict['PFAM_GO_ID'] = pd.concat([pfam_df, go_df], axis=1)

df_dict['labels'] = pd.read_csv('data/prod_activity_identif.csv').iloc[:,1:]

df_merged = pd.merge(df_dict['CDS_motifs'], df_dict['CDS_smCOG'], on = 'BGC_id', how = 'inner')
df_merged = pd.merge(df_merged, df_dict['PFAM_GO_ID'], on = 'BGC_id', how = 'inner')
df_merged = pd.merge(df_merged, df_dict['labels'], on = 'BGC_id', how = 'inner')

df_merged.insert(0, 'cluster_name', df_merged.pop('cluster_name'))
df_merged.insert(1, 'BGC_id', df_merged.pop('BGC_id'))

df_merged.to_csv('data/dataset_merged.csv', index=False)