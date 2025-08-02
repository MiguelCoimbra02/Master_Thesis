#!/usr/bin/python3
#This file has all the functions used to prepare the data to run ConnectTF and also to transform the output of ConnecTF to my use.
# This script is built so that some attributes of the functions can be choses and so, this script can be used for the targeted and the global trimmed networks
# 1st STEP: python3 connecTF.py --function run_before_connecTF --(options) to get both the filter fts and target Network file
# 2nd STEP: Run connecTF with that input and download the target genes list in connecTF download targets option
# 3rd STEP: python3 connecTF.py --function connecTF_network_final --(options) with the connecTF_network, filter, target and loc_at_match files
# 4th STEP: run connecTF with all the 3 inputs: connecTF_network_final, filter and target genes and Additional 3 Edge Features options:
    # ampDAP/Unfiltered
    #DAP/Unfiltered
    #protein-protein interaction
# 5th STEP: python3 connecTF.py --function concatenate_and_transform_connect_tf_parts --(options) to get the final connecTF data table


import pandas as pd
from itertools import product
import os
from pathlib import Path
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Retrieved from another study
blast_file = os.path.join(base_dir, 'connecTF', 'blast_best_hits.txt')
# Script with TF & TR proteins ID extracted from another study
tf_tr_list = os.path.join(base_dir, 'connecTF', 'GCF_002906115.3_Cork_oak_2.0_genomic_longest_isoforms.fasta_ITAK.csv')
network_file = os.path.join(base_dir, 'seidr_output', 'network.csv')
# List obtained in ConnecTF with the TFs available in connecTF to limit my query as much as possible since connecTF does not work with very big queries
TF_available_in_connecTF = os.path.join(base_dir, 'connecTF', 'connecTF_Targeted','latest', 'TF_available_in_connecTF.txt')
output_dir_path = Path(f"{base_dir}/connecTF/connecTF_Targeted/latest")
output_dir_path.mkdir(parents=True, exist_ok=True)
output_dir = os.path.join(base_dir, 'connecTF', 'connecTF_Targeted','latest')


#This function will output 2 files:
    # filter_TFs.txt - File that will be used in the field filter Genes in connecTF
    # LOC_to_Arabidopsis_match.txt - File with 2 cols that matches LOCS to AT genes
def generate_filter_tfs_and_LOC_AT_match(tf_tr_list, blast_file, TF_available_in_connecTF ,output_AT_TF_connecTF, output_AT_LOC_match):
    tf_network_df = pd.read_csv(tf_tr_list, sep=",")
    blast_network_df = pd.read_csv(blast_file, sep="\t")

    #remove TR's
    tf_network_df = tf_network_df[tf_network_df['Type'] == 'TF']  # Keep only rows where Type is 'TF'
    blast_network_df = blast_network_df[blast_network_df['hit'] != '#N/D']  # Remove rows where hit is '#N/D'

    # Only if Protein ID exists in both network_df will be saved
    merged_network_df = tf_network_df.merge(blast_network_df, on="Protein_ID", how="inner")

    # Extract unique Arabidopsis gene names
    arabidopsis_genes = set(merged_network_df["hit"].dropna().unique())
    
    with open(TF_available_in_connecTF, "r") as f:
        available_tf_genes = set(line.strip() for line in f)

    # Keep only genes present in available_AT_tfs_connecTF
    final_genes = arabidopsis_genes.intersection(available_tf_genes)
    
    # Filter Genes To run input to connecTF
    with open(output_AT_TF_connecTF, "w") as f:
        for gene in final_genes:
            f.write(gene + "\n")

     # Create a LOC to Arabidopsis gene mapping (from the merged data)
    loc_mapping = blast_network_df[['Gene_ID', 'hit']].drop_duplicates()

    # Save the LOC to Arabidopsis mapping file
    loc_mapping.to_csv(output_AT_LOC_match, index=False, header=True, sep='\t')

    print(f"Filter Genes file to input in ConnecTF saved: {output_AT_TF_connecTF}")
    print(f"LOC to Arabidopsis mapping saved: {output_AT_LOC_match}.")

    return output_AT_TF_connecTF, output_AT_LOC_match


def generate_Target_Network_ConnecTF(network_file, loc_at_match, output_file):

    # Load the network file
    network_df = pd.read_csv(network_file, delimiter='\t')
    loc_at_df = pd.read_table(loc_at_match, sep='\t')
    loc_at_df.columns = ['Gene_ID', 'hit']
    # Create the new DataFrame for the network
    Target_Network_network_df = pd.DataFrame({
        'source_node': network_df.iloc[:, 0],
        'target_node': network_df.iloc[:, 1],
        'irp': network_df.iloc[:, 2],
    })
    Target_Network_network_df.insert(1, 'interaction', 'interacts_with')
    # Merge network_df with loc_at_df to get homologues for source_node
    merged_df = pd.merge(Target_Network_network_df, loc_at_df, 
                        left_on='source_node', 
                        right_on='Gene_ID', 
                        how='left') \
                .rename(columns={'hit': 'source_homologue'}) \
                .drop(columns=['Gene_ID'])

    # Merge again to get homologues for target_node
    merged_df = pd.merge(merged_df, loc_at_df, 
                        left_on='target_node', 
                        right_on='Gene_ID', 
                        how='left') \
                .rename(columns={'hit': 'target_homologue'}) \
                .drop(columns=['Gene_ID'])

    # Drop rows where either source_homologue or target_homologue is NaN (no homologue found)
    target_network_df = merged_df.dropna(subset=['source_homologue', 'target_homologue'])

    # Replace the source_node and target_node with their respective homologues
    target_network_df.loc[:, 'source_node'] = target_network_df['source_homologue']
    target_network_df.loc[:, 'target_node'] = target_network_df['target_homologue']

    # Drop the temporary columns
    target_network_df = target_network_df.drop(columns=['source_homologue', 'target_homologue'])
    target_network_df.to_csv(output_file, index=False, header=True, sep='\t')
    print(f"Transformation complete. Now run this network in connecTFsaved as {output_file} to fetch the target genes to further reduce the data inputed in connecTF.")


def run_before_connecTF(tf_tr_list = tf_tr_list, blast_file = blast_file , TF_available_in_connecTF = TF_available_in_connecTF, network_file = network_file, output_dir = output_dir):
    '''
    ConnecTF is very bad at handling a lot of data and so this treatment is required to reduce as much as possible the input into connecTF
    In this def im gathering the filter_genes input and forming a network to fetch the target genes input by running connecTF.
    '''
    filter_genes, loc_at_match = generate_filter_tfs_and_LOC_AT_match(tf_tr_list, blast_file, TF_available_in_connecTF, os.path.join(output_dir, 'filter_TFs.txt'), os.path.join(output_dir, 'LOC_to_Arabidopsis_match.txt'))
    generate_Target_Network_ConnecTF(network_file, loc_at_match, os.path.join(output_dir, 'connecTF_network_to_fetch_target_genes.txt'))
    


# ----- After 1st Run of connecTF to fetch target genes -----


targeted_network = os.path.join(base_dir, 'connecTF', 'connecTF_Targeted','latest', 'connecTF_network_to_fetch_target_genes.txt')
target_genes = os.path.join(base_dir, 'connecTF', 'connecTF_Targeted','latest', 'target_genes.txt')
filter_genes = os.path.join(base_dir, 'connecTF', 'connecTF_Targeted','latest', 'filter_TFs.txt')

#Run ConnecTF with TheTarget network to get a list of the targets
def connecTF_network_final(targeted_network = targeted_network, target_genes = target_genes, filter_genes = filter_genes, outfile = os.path.join(output_dir, 'connecTF_network_final.txt')):
    df = pd.read_csv(targeted_network, delimiter='\t')
    tf_arabidopsis_set = set(pd.read_csv(filter_genes, header=None)[0])
    target_set = set(pd.read_csv(target_genes, header=None)[0])

    # Keep rows where a TF is in either 'source_node' or 'target_node' 
    # and the other node is in the target set
    df = df[
        ((df['source_node'].isin(tf_arabidopsis_set)) & (df['target_node'].isin(target_set))) |
        ((df['target_node'].isin(tf_arabidopsis_set)) & (df['source_node'].isin(target_set))) 
    ]
    # Create a list to store the new rows that will be added later
    new_rows = []

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        # Case where target_node is in tf_arabidopsis_set and source_node is in target_set
        if row['target_node'] in tf_arabidopsis_set and row['source_node'] in target_set:
            # Swap the source_node and target_node
            df.at[index, 'source_node'], df.at[index, 'target_node'] = row['target_node'], row['source_node']
        
        # Case where both source_node and target_node are in tf_arabidopsis_set
        elif row['source_node'] in tf_arabidopsis_set and row['target_node'] in tf_arabidopsis_set:
            # Remove the row and create two new ones with swapped values
            df.drop(index, inplace=True)
            new_rows.append({'source_node': row['target_node'], 'target_node': row['source_node']})
            new_rows.append({'source_node': row['source_node'], 'target_node': row['target_node']})

    # Add the new rows to the DataFrame
    df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)
    
    # Save filtered data
    df.to_csv(outfile, index=False, header=True, sep='\t')
    print(f"Transformation complete. File saved as {outfile}.")



#Function used to concat several outputs of connecTF into one big file.
#this was needed because connecTF went down several times when inputing too much data at once

def concat_connec_outputs(output_dir):
    output_parts_dir = os.path.join(output_dir, 'output', 'output_by_parts')
    def is_dir_empty(dir_path):
        return not os.listdir(dir_path)
    
    # Get all CSV files in the directory
    selected_columns = ["gene_id", "TARGET", "EDGE_TYPE", "ADD_EDGES", "gene_name", "TECHNOLOGY/METHOD", "TISSUE/SAMPLE"]
    
    
    
    if is_dir_empty(output_parts_dir):
        # If no partition CSV files are found, it goes look for the full network in the output directory. This happens becasue
        # connecTF sometimes is able to process the full network at once but other times the data must be separated to ease the process.
        # notice that this full network must be called: full_network_connecTF.csv
        
        print("No CSV files found in the directory.")
        df = pd.read_csv(os.path.join(output_dir, 'output', 'output.csv'), usecols=selected_columns)
        print(df.columns.tolist())
        return df
    
    csv_files = [f for f in os.listdir(output_parts_dir) if f.endswith(".csv")]
    # Read the first file with header
    first_file = os.path.join(output_parts_dir, csv_files[0])
    df = pd.read_csv(first_file, usecols=selected_columns)

    for file in csv_files[1:]:
        file_path = os.path.join(output_parts_dir, file)
        try:
            # Attempt to read the CSV file and skip the header row
            temp_df_w_header = pd.read_csv(file_path, usecols=selected_columns) 
            temp_df = temp_df_w_header.iloc[1:].reset_index(drop=True)
            df = pd.concat([df, temp_df], ignore_index=True)
        except Exception as e:
            # Print the file name and the error message
            print(f"Error processing file {file}: {e}")
    
    return df

    #df.to_csv(outfile, index=False)
    #print(f'Final output file from ConnecTF saved as: {outfile}' )

#concat_connec_outputs(connec_out_parts_targeted_dir)


def replace_arabidopsis_with_cork(tf_netw_connecTF, blast_file):
    """
    This function replaces Arabidopsis gene identifiers in a TF network dataset with all corresponding Cork Oak homologs.
    If multiple homologs exist, it generates all possible combinations of interactions.

    Args:
    - tf_netw_connecTF: dataframe object of the TF Arabidopsis network.
    - blast_file: Path to the BLAST results file with Arabidopsis-Cork Oak homolog mappings.
    - output_file: Path to save the transformed TF network CSV.

    Outputs:
    - A CSV file with all possible combinations of homologous interactions.
    """
    
    # Check if input files exist
    '''    if not os.path.exists(tf_netw_connecTF):
        raise FileNotFoundError(f"TF network file not found: {tf_netw_connecTF}")
    if not os.path.exists(blast_file):
        raise FileNotFoundError(f"BLAST file not found: {blast_file}")
    '''    
    # Load TF data
    tf_df = tf_netw_connecTF
    selected_columns = ["gene_id", "TARGET", "EDGE_TYPE", "ADD_EDGES", "gene_name", "TECHNOLOGY/METHOD", "TISSUE/SAMPLE"]
    
    # Check for required columns
    if not all(col in tf_df.columns for col in selected_columns):
        raise ValueError(f"TF network file is missing required columns: {selected_columns}")
    
    tf_df = tf_df[selected_columns]

    # Load BLAST homolog data
    blast_df = pd.read_csv(blast_file, sep="\t", usecols=['Gene_ID', 'hit'])  
    blast_df = blast_df.rename(columns={'Gene_ID': 'cork_gene', 'hit': 'arabidopsis_gene'})

    # Create a mapping where each Arabidopsis gene maps to a list of all its Cork homologs
    arabidopsis_to_cork = blast_df.groupby('arabidopsis_gene')['cork_gene'].apply(list).to_dict()

    # Function to expand a column to all its homologs
    def expand_column(column):
        return column.apply(lambda gene: arabidopsis_to_cork.get(gene, [None]))

    # Expand 'gene_id' and 'TARGET' to all their possible homologs
    tf_df['gene_id_homologs'] = expand_column(tf_df['gene_id'])
    tf_df['target_homologs'] = expand_column(tf_df['TARGET'])

    # Create all possible combinations for rows where both gene_id and TARGET have homologs
    expanded_rows = []
    for _, row in tf_df.iterrows():
        # Get all possible homolog combinations for this row
        homolog_combinations = product(row['gene_id_homologs'], row['target_homologs'])
        for gene_id_homolog, target_homolog in homolog_combinations:
            if gene_id_homolog and target_homolog:  # Ensure neither is None
                new_row = row.copy()
                new_row['gene_id'] = gene_id_homolog
                new_row['TARGET'] = target_homolog
                expanded_rows.append(new_row)

    # Create the final DataFrame with expanded rows
    expanded_df = pd.DataFrame(expanded_rows, columns=selected_columns)

    # Save the transformed network

    return expanded_df
    #expanded_df.to_csv(output_file, index=False, header=True, sep="\t")
    #print(f"Transformation complete. File saved as {output_file}.")

#replace_arabidopsis_with_cork(connec_out_targeted, blast_file)

def clean_connec_output_table(connecTF_table=os.path.join(output_dir, 'output', 'QS_output.csv'), output_file=os.path.join(output_dir, 'output', 'QS_output_final.csv')):
        # Read the table
    connecTF_table = connecTF_table

    # Group by 'gene_id' and 'TARGET' and concatenate 'EDGE_TYPE' values with '/' as separator
    concatenated = connecTF_table.groupby(['gene_id', 'TARGET'], as_index=False).agg({
        'EDGE_TYPE': lambda x: ' / '.join(x.unique())  # Concatenate EDGE_TYPE values with '/' separator
    })

    # Drop the EDGE_TYPE column from connecTF_table before the merge
    connecTF_table_cleaned = connecTF_table.drop(columns=['EDGE_TYPE'])

    # Perform the merge, including only the EDGE_TYPE from concatenated
    result = pd.merge(
        connecTF_table_cleaned,
        concatenated[['gene_id', 'TARGET', 'EDGE_TYPE']],  # Select only the necessary columns
        on=['gene_id', 'TARGET'],
        how='left'
    )

    # Drop duplicates to only keep the first occurrence
    result = result.drop_duplicates(subset=['gene_id', 'TARGET'])

    # Print the resulting DataFrame
    print(result.sample(n=50))
    result.to_csv(output_file, index=False, header=True, sep="\t")

#clean_connec_output_table()

def concatenate_and_transform_connect_tf_parts(output_parts_dir = os.path.join(output_dir, 'output', 'output_by_parts'), blast_file = blast_file, output_file=os.path.join(output_dir, 'output', 'QS_output_final.csv')):
    concatenated_output = concat_connec_outputs(output_parts_dir)
    replaced_df = replace_arabidopsis_with_cork(concatenated_output, blast_file)
    clean_connec_output_table(replaced_df, output_file)





import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Utility script to run ConnecTF processing steps.",
        epilog="""
                Examples:

                Run preprocessing:
                    python run_connect_tf.py --function run_before_connecTF --blast_file path/to/blast.tsv --network path/to/network.csv --output_dir results/

                Run final network creation:
                    python run_connect_tf.py --function connecTF_network_final --network path/to/network.csv --target_genes genes.txt --filter_genes filters.txt

                Concatenate ConnecTF output parts:
                    python run_connect_tf.py --function concatenate_and_transform_connect_tf_parts --output_dir results/ --blast_file path/to/blast.tsv

                Use --help with any argument for more info.
                """
    )

    parser.add_argument('--function', required=True, choices=[
        'run_before_connecTF',
        'connecTF_network_final',
        'concatenate_and_transform_connect_tf_parts'
    ], help='Function to run.')

    # Optional arguments depending on the function
    parser.add_argument('--output_dir', type=str, default=f'{os.path.join(base_dir, 'connecTF', 'connecTF_Targeted','latest')}', help='Path to output directory.')
    parser.add_argument('--blast_file', type=str, default=f'{os.path.join(base_dir, 'connecTF', 'blast_best_hits.txt')}', help='Path to BLAST file.')
    parser.add_argument('--target_genes', type=str, default=f'{os.path.join(base_dir, 'connecTF', 'connecTF_Targeted','latest', 'target_genes.txt')}',help='Path to target genes file.')
    parser.add_argument('--filter_genes', type=str, default=f'{os.path.join(base_dir, 'connecTF', 'connecTF_Targeted','latest', 'filter_TFs.txt')}',help='Path to filter genes file.')
    parser.add_argument('--network', type=str, default=f'{os.path.join(base_dir, 'seidr_output', 'network.csv')}',help='Path to network file.')
    

    args = parser.parse_args()

    if args.function == "run_before_connecTF":
        run_before_connecTF(
            tf_tr_list,
            blast_file,
            TF_available_in_connecTF,
            args.network,
            args.output_dir
        )

    elif args.function == "connecTF_network_final":
        connecTF_network_final(args.network, args.target_genes, args.filter_genes, args.output_dir)

    elif args.function == "concatenate_and_transform_connect_tf_parts":
        concatenate_and_transform_connect_tf_parts(args.output_dir, args.blast_file, f'{args.output_dir}/output/QS_output_final.csv')

    else:
        parser.print_help()


if __name__ == "__main__":
    main()