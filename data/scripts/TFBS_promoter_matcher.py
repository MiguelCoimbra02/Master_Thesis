#!/usr/bin/python3
import pandas as pd
import argparse
"""
This script will be runned in fimo.sh
This script filters a co-expression network based on transcription factor (TF) IDs, extracts matching promoter sequences, and outputs them in a FASTA format.
needs the following files:
1) co expression network file w/ source/target/gene cols
2) promotor sequences file w/ Sequence/Scaffold/Gene
3) TF in quercus format
"""
def TFBS_promoter_matcher(network, promoter_seq, tf, out_dir):
    
    network = pd.read_csv(network, sep="\t")
    query_promoter_sequences = pd.read_csv(promoter_seq, sep=",")
    target_list=[]
    filtered_network = network[network["Source"].str.contains(tf) | network["Target"].str.contains(tf)]

    for gene in filtered_network["Source"]:
        if gene != tf:
            target_list.append(gene)
    
    for gene in filtered_network["Target"]:
        if gene != tf:
            target_list.append(gene)

    queried_rows = query_promoter_sequences.loc[query_promoter_sequences["Gene"].isin(target_list)] 
    fasta_output = open(out_dir + "/TFBS_TF_" + tf + "_fasta.fa", "w")
    for row_index in range(0, len(queried_rows)):
        S = queried_rows["Scaffold"].values[row_index]
        G = queried_rows["Gene"].values[row_index]
        seq = queried_rows["Sequence"].values[row_index]
        one_entry = ">" + G + " " + S + "\n" + seq + "\n"
        one_entry = one_entry.upper()
        fasta_output.write(one_entry)
    fasta_output.close()

def main():
    parser = argparse.ArgumentParser(description="Filters a co-expression network based on TF IDs, extracts matching promoter sequences, and outputs them in FASTA format.")
    parser.add_argument('--network', required=True, metavar="FILE", help="Co-expression network file (tab-separated)")
    parser.add_argument('--promoter_seq', required=True, metavar="FILE", help="Promoter sequences file (comma-separated)")
    parser.add_argument('--tf', required=True, metavar="STR", help="TF ID in Quercus format")
    parser.add_argument('--out_dir', type=str, default='/home/mlcoimbraj/data/cis_elements/TF_motifs/fasta', help="Output directory")

    args = parser.parse_args()

    TFBS_promoter_matcher(args.network, args.promoter_seq, args.tf, args.out_dir)

if __name__ == "__main__":
    main()


