#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This Script loads an cor_scores.txt with the correlation scores and outputs a cor_scores_mean.txt in a ready format to import into Cytoscape
"""

import pandas as pd
import argparse

def cor_scores_mean(cor_scores, outname):
    df = pd.read_table(cor_scores, sep='\t')  # Loads the correlation scores file into a pandas DataFrame.
    print(df.head())
    print('--------------------')
    print(df.columns)
    df['Gene_Pair'] = df['Source'] + ' (interacts with) ' + df['Target']  # Creates a 'Gene_Pair' column.
    df_cor = df.assign(cor_mean=lambda x: df.iloc[:, 2:5].mean(axis=1))  # Computes the mean of the correlation scores across columns 2-4.
    df_cor.to_csv(outname)  # Saves the result to the output file.

def main():
    parser=argparse.ArgumentParser(description='Add a mean column to the cor_scores table')
    parser.add_argument('--cor_scores', required=True, metavar='tab', 
                        help='cor_scores for which it will be calculated the mean by row')
    parser.add_argument('--out', metavar='STR', type=str, default="cor_scores_mean.txt",
                        help="Output file name")
    
    args = parser.parse_args()
    
    cor_scores_mean(args.cor_scores, args.out)

    
main()