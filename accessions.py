import pandas as pd

acc_antifungi = pd.read_csv('data/antifungi/all_features_CDS_motifs.csv')['acc']
acc_nonantifungi = pd.read_csv('data/antifungi/all_features_CDS_motifs.csv')['acc']

accessions = pd.concat([acc_antifungi, acc_nonantifungi], axis = 0)

print(accessions.unique())