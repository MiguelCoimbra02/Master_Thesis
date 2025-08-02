#!/usr/bin/python3
# >>> python3 extract_TF_TR_from_global_network.py
# This script extracts the TF-TF, TF-TR, TR-TF and TR-TR interactions from the global network.
# This is necessary since when running seidr in Targeted mode, these interactions are not included in the network.
import pandas as pd
import os
from pathlib import Path

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
global_network_dir = Path(os.path.join(base_dir, 'seidr_output_0_05'))
blast_file = Path(os.path.join(base_dir, 'connecTF', 'blast_best_hits.txt'))
connecTF_files_dir = Path(os.path.join(base_dir, 'connecTF', 'connecTF_Targeted', 'latest'))
output_dir = global_network_dir

def extract_TF_TR_from_global_network(global_network = f'{global_network_dir}/network.csv', tf_tr_file = f'{connecTF_files_dir}/tf_and_tr_file.txt', output_file = f'{output_dir}/trimmed_network.csv'):
    # Load data
    global_network_df = pd.read_csv(global_network, sep='\t')

    # Set with genes in the network
    with open(tf_tr_file, 'r') as genes_file:
        genes_of_interest = set(genes_file.read().splitlines())
    
    # Filter rows where both Source and Target are in genes_of_interest
    filtered_df = global_network_df[
        global_network_df['Source'].isin(genes_of_interest) &
        global_network_df['Target'].isin(genes_of_interest) &
        (global_network_df['irp_score'] >= 0.15) # minimum irp_score present in my targeted dataset (rounded up)
    ]

    # Save the trimmed network
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"Trimmed network with only TF-TF, TF-TR, TR-TF and TR-TR interactions saved to {output_file} with {len(filtered_df)} rows.")

extract_TF_TR_from_global_network()



