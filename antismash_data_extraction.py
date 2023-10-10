import pandas as pd
import re
import os
from Bio import SeqIO



######################################################################
### Cria um arquivo para adicionar os nomes e numeros dos clusters ###
######################################################################

clusters_names_csv = open('clustersname.csv','w')

for dirpath, dirnames, files in os.walk('./'):
    for file_name in files:
        file_path = os.path.join(dirpath, file_name)
        
        ### busca arquivos dentro da pasta /knownclusterblast, output do antiSMASH
        if dirpath == './knownclusterblast':
            
            ### Obtem o numero de acesso da cepa que foi analisada pelo antiSMASH
            
            acc_clust = file_name.split('.')[0]
            
            ### extrai o numero do cluster
            
            str_fn = str(file_name).replace('.txt','').replace('c','').split('_')[-1]
            
            ### abre o arquivo de texto e le as linhas
            
            file = open(file_path,'r')            
            lines = file.readlines()
            cluster_name = ''
            
            ### Corre pelas linhas do arquivo com enumerate para conseguir manipular indexes das linhas
            
            for index, lin in enumerate(lines):
                if lin.startswith('Significant hits'):
                    if index > 0:
                        previous = lines[index - 1]
                    if index < (len(lines) - 1):
                        next_ = lines[index + 1]
                        
                        ### obtem o nome (caso tenha) do cluster, na linha seguinte a Sginificant Hits
                        
                        if next_.startswith('1.'):
                            cluster_name = next_.replace('\t','_').split('_')[-1]
                        else:
                            cluster_name = 'region'
            clusters_names_csv.write(str_fn+','+acc_clust+',')
            if cluster_name != None:
                clusters_names_csv.write(str(cluster_name)+'\n')
        else:
            pass
clusters_names_csv.close()


##############################################
### Tratamento do CSV criado anteriormente ###
### pulando linhas e espaços em branco     ###
##############################################

clust_name = pd.read_csv('clustersname.csv', 
                         skip_blank_lines=True, 
                         header=None)

### Renomeia as colunas

clust_name.rename(columns={0:'clust_num',1:'n_acesso',2:'nome'}, inplace = True)


### Coloca em ordem com base no numero do cluster
df2 = clust_name.sort_values(by=['n_acesso','clust_num'])

### adiciona a string "region" nos clusters onde nao foi possivel obter o nome
df2.fillna('region', inplace = True)

### Trata os indexes
df2.reset_index(drop=True, inplace=True)
df2.index = df2.index + 1
print(df2)
        
clusters_names_csv.close()


################################################
### Leitura dos arquivos output do antiSMASH ###
################################################

os.mkdir('data_frame/')

for dirpath, dirnames, files in os.walk('./'):
    
    for dirname in dirnames:
        pastas = dirname
    for file_name in files:
        
        ### Seleciona os arquivos dos clusters output do antiSMASH
        
        if file_name.startswith('NZ' or 'CP') and file_name.endswith('.gbk'):
            
            cds_motifs_list = []
            cds_motifs_counts = {}
            
            CDS_smCOG_list = []
            CDS_smCOG_count = {}
            
            PFAM_desc_list = []
            PFAM_desc_count = {}  
            
            PFAM_GO_list = []
            PFAM_GO_count = {}
            
            PFAM_ID_list = []
            PFAM_ID_count = {}
            
            ### Obtem o numero do cluster dentro do resultado do antiSMASH
            
            reg_num0 = file_name.split('.')[-2][-2:]
            
            ### Obtem o numero de acesso da cepa que foi analisada pelo antiSMASH
            
            acc_clust = file_name.split('.')[0]
                        
            ### Abre o arquivo
            
            record = SeqIO.read(file_name, "genbank")
            
            
            ### Corre sobre as features do arquivo gbk
            for features in record.features:
                
                ### Trata e extrai as informações de features do tipo CDS_motifs
                ### Usa os qualifiers 'note' e 'aSTool'
                if features.type == 'CDS_motif':
                    if features.qualifiers.get('note') != None:
                        cdsmotif_note = features.qualifiers.get('note')[0]
                        cdsmotif_note_f = cdsmotif_note.split(':')[1].strip()
                        if cdsmotif_note_f not in cds_motifs_list:
                            cds_motifs_list.append(cdsmotif_note_f)
                            cds_motifs_counts[cdsmotif_note_f] = 0
                        cds_motifs_counts[cdsmotif_note_f] += 1
                        
                    elif features.qualifiers.get('aSTool') != None:
                        #print(features.qualifiers.get('aSTool'))
                        labels_motifs = str(features.qualifiers.get('label')).strip('[]"\'"')
                        if labels_motifs not in cds_motifs_list:
                            cds_motifs_list.append(labels_motifs)
                            #print(labels_motifs)
                            cds_motifs_counts[labels_motifs] = 0
                        cds_motifs_counts[labels_motifs] += 1                        
                        
                    else:
                        cdsmotif_note = None 
                        cdsmotif_note_f = None
                
                if features.type == 'CDS':
                    gen_func = features.qualifiers.get('gene_functions')
                    if gen_func != None:
                        for items in gen_func:                        
                            if 'smcogs' in items:
                                SMCOG_type0 = re.search("SMCOG[0-9]*:.*[(]", str(items))
                                SMCOG_type1 = SMCOG_type0.group().replace("(","").strip()
                                if SMCOG_type1 not in CDS_smCOG_list:
                                    CDS_smCOG_list.append(SMCOG_type1)
                                    CDS_smCOG_count[SMCOG_type1] = 0
                                CDS_smCOG_count[SMCOG_type1] += 1
                
                if features.type == 'PFAM_domain':
                    
                    PFAM_desc = str(features.qualifiers.get('description'))
                    PFAM_desc2 = PFAM_desc.strip('[]"\'"')
                    #print(PFAM_desc2)
                    if PFAM_desc2 not in PFAM_desc_list:
                        PFAM_desc_list.append(PFAM_desc2)
                        PFAM_desc_count[PFAM_desc2] = 0
                    PFAM_desc_count[PFAM_desc2] += 1
                    
                    GO_terms = features.qualifiers.get('gene_ontologies')
                    if GO_terms != None:
                        db_xref = features.qualifiers.get('db_xref')                        
                        for refs in db_xref:
                            if refs.startswith('GO:'):
                                GO = refs
                                if GO not in PFAM_GO_list:
                                    PFAM_GO_list.append(GO)
                                    PFAM_GO_count[GO] = 0
                                PFAM_GO_count[GO] += 1
                                
                    pfam_id = features.qualifiers.get('db_xref')[0]
                    if pfam_id not in PFAM_ID_list:
                        PFAM_ID_list.append(pfam_id)
                        PFAM_ID_count[pfam_id] = 0
                    PFAM_ID_count[pfam_id] += 1
                    #print(pfam_id)
                           
            
            ### Cria 1 DataFrame pra cada cluster contendo a contagem de cada feature extraído
            
            df1 = pd.DataFrame([cds_motifs_counts], index=[1])
            
            ### Extrai nome do cluster baseado no seu numero de acordo com a DataFrame dos clusters_names
            
            cluster_name_ = df2.loc[(df2["clust_num"] == int(reg_num0)) 
                                     & (df2['n_acesso'] == acc_clust) , "nome"]
            
            ### transforma a Series obtida em DataFrame
            
            cf = cluster_name_.to_frame()
            
            ### Adiciona a coluna 'nome' com o nome do cluster obtido anteriormente e a coluna 'acc' com numero de acesso da cepa
            
            df1['nome'] = cf.iloc[0,0]
            df1['acc'] = acc_clust
            
            #print(df1)
            
            ### Salva o arquivo numa pasta de CSV de cada cluster
            df1.to_csv('data_frame/'+acc_clust+'tempfile_motifs_cluster'+reg_num0+'.csv')
            
            
            cds_df = pd.DataFrame([CDS_smCOG_count], index=[1])
            cluster_name_1 = df2.loc[(df2["clust_num"] == int(reg_num0)) 
                                     & (df2['n_acesso'] == acc_clust) , "nome"]
            cf1 = cluster_name_1.to_frame()
            cds_df['nome'] = cf1.iloc[0,0] 
            cds_df['acc'] = acc_clust
            cds_df.to_csv('data_frame/'+acc_clust+'tempfile_smCOG'+reg_num0+'.csv')
            
            PFAM_df = pd.DataFrame([PFAM_desc_count], index=[1])
            cluster_name_2 = df2.loc[(df2["clust_num"] == int(reg_num0)) 
                                     & (df2['n_acesso'] == acc_clust) , "nome"]
            cf2 = cluster_name_2.to_frame()
            PFAM_df['nome'] = cf2.iloc[0,0]
            PFAM_df['acc'] = acc_clust
            PFAM_df.to_csv('data_frame/'+acc_clust+'tempfile_PFAMdesc'+reg_num0+'.csv')
            
            PFAM_GO = pd.DataFrame([PFAM_GO_count], index=[1])
            cluster_name_3 = df2.loc[(df2["clust_num"] == int(reg_num0)) 
                                     & (df2['n_acesso'] == acc_clust) , "nome"]
            cf3 = cluster_name_3.to_frame()
            PFAM_GO['nome'] = cf3.iloc[0,0]
            PFAM_GO['acc'] = acc_clust
            PFAM_GO.to_csv('data_frame/'+acc_clust+'tempfile_PFAM_GO'+reg_num0+'.csv')
            
   
            PFAM_ID = pd.DataFrame([PFAM_ID_count], index=[1])
            cluster_name_4 = df2.loc[(df2["clust_num"] == int(reg_num0)) 
                                     & (df2['n_acesso'] == acc_clust) , "nome"]
            cf4 = cluster_name_4.to_frame()
            PFAM_ID['nome'] = cf4.iloc[0,0]
            PFAM_ID['acc'] = acc_clust
            PFAM_ID.to_csv('data_frame/'+acc_clust+'tempfile_PFAM_ID'+reg_num0+'.csv')
            
####################################################            
### Concatena todos os CSV criados eu um arquivo ###
####################################################

txt0 = 'motifs_cluster'

files_paths0 = []

for dirpath, dirnames, files in os.walk('./'):
    for file_name in files:
        file_path0 = os.path.join(dirpath, file_name)
        if dirpath == './data_frame' and txt0 in file_name:
            
            files_paths0.append(file_path0)

df5 = pd.concat(map(pd.read_csv, files_paths0), ignore_index=True)
df5 = df5.fillna(0)
df5 = df5.sort_values(by='acc')
df5.to_csv('data_frame/all_features_CDS_motifs.csv', header=True)


txt1 = 'smCOG'

files_paths1 = []

for dirpath, dirnames, files in os.walk('./'):
    for file_name in files:
        file_path1 = os.path.join(dirpath, file_name)
        if dirpath == './data_frame' and txt1 in file_name:
            
            files_paths1.append(file_path1)
            
df6 = pd.concat(map(pd.read_csv, files_paths1), ignore_index=True)
df6 = df6.fillna(0)
df6 = df6.sort_values(by='acc')
df6.to_csv('data_frame/all_features_CDS_smCOG.csv', header=True)


txt2 = 'PFAMdesc'

files_paths2 = []

for dirpath, dirnames, files in os.walk('./'):
    for file_name in files:
        file_path2 = os.path.join(dirpath, file_name)
        if dirpath == './data_frame' and txt2 in file_name:
            
            files_paths2.append(file_path2)
            
df7 = pd.concat(map(pd.read_csv, files_paths2), ignore_index=True)
df7 = df7.fillna(0)
df7 = df7.sort_values(by='acc')
df7.to_csv('data_frame/all_features_PFAM_desc.csv', header=True)


txt3 = 'PFAM_GO'

files_paths3 = []

for dirpath, dirnames, files in os.walk('./'):
    for file_name in files:
        file_path3 = os.path.join(dirpath, file_name)
        if dirpath == './data_frame' and txt3 in file_name:
            
            files_paths3.append(file_path3)
            
df8 = pd.concat(map(pd.read_csv, files_paths3), ignore_index=True)
df8 = df8.fillna(0)
df8 = df8.sort_values(by='acc')
df8.to_csv('data_frame/all_features_PFAM_GO.csv', header=True)


txt4 = 'PFAM_ID'

files_paths4 = []

for dirpath, dirnames, files in os.walk('./'):
    for file_name in files:
        file_path4 = os.path.join(dirpath, file_name)
        if dirpath == './data_frame' and txt4 in file_name:
            #print(file_name)
            files_paths4.append(file_path4)
            
df9 = pd.concat(map(pd.read_csv, files_paths4), ignore_index=True)
df9 = df9.fillna(0)
df9 = df9.sort_values(by='acc')
df9.to_csv('data_frame/all_features_PFAM_ID.csv', header=True)

for dirpath, dirnames, files in os.walk('./'):
    for file_name in files:
        file_path__ = os.path.join(dirpath, file_name)
        if re.search('tempfile', file_name):
            os.remove(file_path__)
