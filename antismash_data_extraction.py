import pandas as pd
import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


######################################################################
### Cria um arquivo para adicionar os nomes e numeros dos clusters ###
######################################################################

clusters_names_csv = open('clustersname.csv','w')

for dirpath, dirnames, files in os.walk('./'):
    for file_name in files:
        file_path = os.path.join(dirpath, file_name)
        
        ### busca arquivos dentro da pasta /knownclusterblast, output do antiSMASH
        if dirpath == './knownclusterblast':
            
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
            clusters_names_csv.write(str_fn+',')
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

clust_name.rename(columns={0:'clust_num',1:'nome'}, inplace = True)


### Coloca em ordem com base no numero do cluster
df2 = clust_name.sort_values(by='clust_num')

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

for dirpath, dirnames, files in os.walk('./'):
    for dirname in dirnames:
        pastas = dirname
    for file_name in files:
        
        ### Seleciona os arquivos dos clusters output do antiSMASH
        
        if file_name.startswith('NZ_CP059318.1.region0') and file_name.endswith('.gbk'):
            
            cds_motifs_list = []
            cds_motifs_counts = {}
            
            
            ### Obtem o numero do clusters dentro do resultado do antiSMASH
            
            reg_num0 = file_name.split('.')[-2][-2:]
            
            ### Obtem o numero de acesso da cepa que foi analisada pelo antiSMASH
            
            acc_clust = file_name.split('.')[0]
                        
            ### Abre o arquivo
            
            record = SeqIO.read(file_name, "genbank")
            #print('\n'+file_name)
            
            
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
            
            ### Cria 1 DataFrame pra cada cluster contendo a contagem de cada feature extraído
            
            df1 = pd.DataFrame([cds_motifs_counts], index=[1])
            
            ### Extrai nome do cluster baseado no seu numero de acordo com a DataFrame dos clusters_names
            
            cluster_name_ = df2.loc[df2["clust_num"] == int(reg_num0), "nome"]
            
            ### transforma a Series obtida em DataFrame
            
            cf = cluster_name_.to_frame()
            
            ### Adiciona a coluna 'nome' com o nome do cluster obtido anteriormente e a coluna 'acc' com numero de acesso da cepa
            
            df1['nome'] = cf.iloc[0,0]
            df1['acc'] = acc_clust
            
            ### Salva o arquivo numa pasta de CSV de cada cluster
            df1.to_csv('data_frame/tempfile'+reg_num0+'.csv')
            
####################################################            
### Concatena todos os CSV criados eu um arquivo ###
####################################################

files_paths = []

for dirpath, dirnames, files in os.walk('./'):
    #print(f'Found directory: {dirpath}')
    for file_name in files:
        file_path = os.path.join(dirpath, file_name)
        if dirpath == './data_frame' and file_name.startswith('tempfile'):
            
            files_paths.append(file_path)
            
df5 = pd.concat(map(pd.read_csv, files_paths), ignore_index=True)
df5 = df5.fillna(0)
df5.to_csv('data_frame/all_features_CDS_motifs.csv', header=True)
