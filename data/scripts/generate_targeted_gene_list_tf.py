#!/usr/bin/python3
# This script creates files with the TF and TR Quercus genes, 1 per line, to run Seidr in targeted mode.
# It uses 3 files:
# 1 - List of proteins of identified Quercus TFs
# 2 - Matches between proteins and Quercus genes
# 3 - genes.txt file to be used as input in Seidr

import pandas as pd
import os
from pathlib import Path

# Base directory setup
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
tf_protein_list = Path(os.path.join(base_dir, 'connecTF', 'GCF_002906115.3_Cork_oak_2.0_genomic_longest_isoforms.fasta_ITAK.csv'))
blast_file = Path(os.path.join(base_dir, 'connecTF', 'blast_best_hits.txt'))
genes_seidr_file = Path(os.path.join(base_dir, 'seidr_output', 'genes.txt'))
out_dir = Path(os.path.join(base_dir, 'connecTF', 'connecTF_Targeted', 'latest'))

# Load data
tf_protein_df = pd.read_csv(tf_protein_list, sep=',')
blast_file_df = pd.read_csv(blast_file, sep='\t')

# Set with genes in the network
with open(genes_seidr_file, 'r') as genes_file:
    genes_on_network_set = set(genes_file.read().splitlines())

# Merge and get all TF/TR genes
tf_tr_merged = pd.merge(tf_protein_df, blast_file_df, on='Protein_ID')
loc_tf_tr_set_all = set(tf_tr_merged['Gene_ID'].tolist())

# Filter to keep only genes present in the network
loc_tf_tr_set = {gene for gene in loc_tf_tr_set_all if gene in genes_on_network_set}

# Write TF and TR combined file
out_dir.mkdir(parents=True, exist_ok=True)
with open(out_dir / 'tf_and_tr_file.txt', 'w') as f:
    for gene in loc_tf_tr_set:
        f.write(f"{gene}\n")
print(f'File with TF and TR that exist in the network written to {out_dir} directory')

# TF File LOCs
tf_merged = tf_tr_merged[tf_tr_merged['Type'] == 'TF']
loc_tf_set_all = set(tf_merged['Gene_ID'].tolist())
loc_tf_set = {gene for gene in loc_tf_set_all if gene in genes_on_network_set}

with open(out_dir / 'tf_file.txt', 'w') as f:
    for gene in loc_tf_set:
        f.write(f"{gene}\n")
print(f'File with only TFs that exist in the network written to {out_dir} directory')

# TR File LOCs
tr_merged = tf_tr_merged[tf_tr_merged['Type'] == 'TR']
loc_tr_set_all = set(tr_merged['Gene_ID'].tolist())
loc_tr_set = {gene for gene in loc_tr_set_all if gene in genes_on_network_set}

with open(out_dir / 'tr_file.txt', 'w') as f:
    for gene in loc_tr_set:
        f.write(f"{gene}\n")
print(f'File with only TRs that exist in the network written to {out_dir} directory')
