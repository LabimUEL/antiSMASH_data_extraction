import pandas as pd
import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

for dirpath, dirnames, files in os.walk('./'):
    for dirname in dirnames:
        pastas = dirname
    for file_name in files:
        
        if file_name.startswith('NZ_CP059318.1.region0') and file_name.endswith('.gbk'):
            
            cds_motifs_list = []
            cds_motifs_counts = {}
            
            record = SeqIO.read(file_name, "genbank")
            print('\n'+file_name)
            for features in record.features:
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
                        
                    CDS_motif = open('features_matrices/CDS_motifs.csv', 'w')
                    
                   
                    

            for i in cds_motifs_counts:
                for j in cds_motifs_list:
                    if j in cds_motifs_counts:
                        CDS_motif.write(str(cds_motifs_counts[i]) + ",")
                    else:
                        CDS_motif.write("0,")
                CDS_motif.write("\n")        
            print(cds_motifs_counts)
            print(len(cds_motifs_list))
        CDS_motif.close()