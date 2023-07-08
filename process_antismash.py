import pandas as pd

df_acc = pd.read_csv('data/antifungi_acc.csv')

datasets = {
    'df_cds_motifs': pd.read_csv('data/antismash/ALL_features_CDS_motifs.csv'),
    'df_cds_smcog': pd.read_csv('data/antismash/ALL_features_CDS_smCOG.csv'),
    'df_pfam_desc': pd.read_csv('data/antismash/ALL_features_PFAM_desc.csv'),
    'df_pfam_go': pd.read_csv('data/antismash/ALL_features_PFAM_GO.csv'),
    'df_pfam_id': pd.read_csv('data/antismash/ALL_features_PFAM_ID.csv')
}

for key, df in datasets.items():
    df['label'] = df['acc'].str.contains('|'.join(df_acc['acc'].astype(str)), na=False).astype(int)
    
    # Change the position of the 'acc' column to the last column
    cols = df.columns.tolist()
    cols.remove('acc')
    cols.append('acc')
    cols.remove('nome')
    cols.append('nome')
    df = df.reindex(columns=cols)
    # Update the DataFrame in the dictionary
    datasets[key] = df.iloc[:, 3:]

# Output the updated DataFrames
for key, df in datasets.items():
    df.to_csv(f"data/processed/{key}.csv", index=False)