import pandas as pd

df1 = pd.read_csv('data/all_features_PFAM_ID.csv').iloc[:,2:]
df2 = pd.read_csv('data/prod_activity_identif.csv').iloc[:,1:]

df3 = pd.merge(df1, df2, how = 'inner', on = 'BGC_id')

df3.insert(0, 'BGC_id', df3.pop('BGC_id'))

df3.to_csv('data/dataset_merged.csv', index=False)