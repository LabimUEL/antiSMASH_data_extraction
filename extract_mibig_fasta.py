from Bio import SeqIO, SearchIO
from numpy import loadtxt
import pandas as pd
import subprocess

file = open('data/BGCs_com_atividade.txt', 'rb')
bgcs = list(loadtxt(file, dtype = str, delimiter = ","))
bgcs.sort()

print(bgcs)

def extract_activity_BGCs():
    with open('data/mibig_prot_seqs_2.0.fasta') as mibig_fasta:
        for seq_record in SeqIO.parse(mibig_fasta, 'fasta'):
            if any(ext in seq_record.description for ext in bgcs):
                print(seq_record.description)

                with open('data/mibig_seqs_activity.fasta', 'a') as activity_fasta:
                    activity_fasta.write('>' + seq_record.description + '\n')
                    activity_fasta.write(str(seq_record.seq) + '\n')
                    activity_fasta.close()

        mibig_fasta.close()

def hmmer_table():
    subprocess.run(['sh', 'hmm_extract.sh'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    
    df = pd.read_csv('data/pfam.tsv', delimiter = '\t', names = ['target name', 'accession', 'query name', 'description of target'])

    df.to_csv('data/pfam.tsv', sep ='\t', index = False)

    print(df)