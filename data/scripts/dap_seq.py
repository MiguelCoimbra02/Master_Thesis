#!/usr/bin/python3
#This file will retrieve the targets for the genes specified in the name of each file and turn them into a df
#with collumn gene_id and TARGET so that in a later script this can be merged with the network df
#Requires a file with the gene name as title and its targets one per line
#RUN: python3 dap_seq.py LOC1.txt LOC2.txt ... 
#I runned it like this: python3 dap_seq.py LOC111997151.txt LOC112008346.txt LOC112030452.txt

import pandas as pd
import os
import argparse
from pathlib import Path

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
dap_seq_dir = Path(os.path.join(base_dir, 'dap_seq')) 
os.makedirs(dap_seq_dir, exist_ok=True) #exist_ok to avoid errors if the dir already exists


parser = argparse.ArgumentParser(description="Read multiple text files and turn them into a dataframe with columns gene_id and TARGET.")
parser.add_argument('files', type=str, nargs='+', help='Names of the text files where the name is the Gene and the content its targets, one per line')
args = parser.parse_args()

#will create the df from the dictionarys (rows) here stored
data = []

for file in args.files:
    file_path = dap_seq_dir / file
    with open(file_path, 'r') as gene_file:
        targets = gene_file.read().splitlines()
        gene_id = file_path.stem  # Get the file name without extension as gene_id

        for target in targets:
            data.append({'gene_id': gene_id, 'TARGET': target})
            # reverse interaction to later merge with network table
            data.append({'gene_id': target, 'TARGET': gene_id})

dap_seq_df = pd.DataFrame(data).drop_duplicates()
dap_seq_df.to_csv(f'{dap_seq_dir}/dap_seq.txt', sep='\t', index=False)