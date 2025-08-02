#!/usr/bin/python3
# runned this for targeted:python generate_edges_nodes_tables.py --function generate_edges_table --network ../seidr_output/network.csv --out_dir ../data/data_targeted/ --connecTF_table ../connecTF/connecTF_Targeted/latest/output_targeted/QS_output_final.csv  --output_edges_name edge_full.csv
# runned this for TF-TF interactions: python generate_edges_nodes_tables.py --function generate_edges_table --network ../seidr_output_0_05/trimmed_network.csv --out_dir ../data/ --connecTF_table ../connecTF/output/QS_output_final.csv  --output_edges_name edges_TF_TR_final.csv

import pandas as pd
import os
from pathlib import Path


'''base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
network_dir = Path(os.path.join(base_dir, 'seidr_output'))
connecTF_files_dir = Path(os.path.join(base_dir, 'connecTF','connecTF_Targeted', 'latest'))
connecTF_dir = Path(os.path.join(connecTF_files_dir, 'output_targeted'))
#node_specific 
loc_to_at_file = Path(os.path.join(connecTF_files_dir, 'LOC_to_Arabidopsis_match.txt'))
tf_file = Path(os.path.join(connecTF_files_dir, 'tf_file.txt'))
tr_file = Path(os.path.join(connecTF_files_dir, 'tr_file.txt'))
gene_annotation_file = Path(os.path.join(base_dir, 'genome', 'genome_annotation', 'gene_functional_annotation.csv'))

#edge_specific
cis_dir = Path(os.path.join(base_dir, 'cis_elements', 'fimo_out'))
dap_seq_dir = Path(os.path.join(base_dir, 'dap_seq')) #need to create a dir to keep the files with dap_seq targets for the 3 genes
out_dir = Path(os.path.join(base_dir, 'data', 'data_targeted', 'latest'))
os.makedirs(out_dir, exist_ok=True) #exist_ok to avoid errors if the dir already exists
'''
#--------------- utils ---------------
    
def determine_tf_rank(row):
    tf_rank = [0]  # co-exp -> 0
    
    if row['connecTF'] == 1:
        tf_rank.append(1)  # connecTF -> 1
    if row['cis_elements'] == 1:
        tf_rank.append(2)  # cis_elements -> 2
    if row['dap_seq'] == 1:
        tf_rank.append(3)  # dap_seq -> 3
    return (tf_rank)


def determine_direction(row):
    directed = 0
    if len(row['tf_rank']) > 1:
        directed = 1
    return directed
    

def generate_nodes_table(network_table, tf_file, tr_file, gene_annotation_file, loc_to_at_file, output_directory, output_file_name):
    '''
        name VARCHAR(255) PRIMARY KEY,
         Arabidopsis_gene VARCHAR(255),
         isTF INT,
         isTR INT,
         gene_annotation VARCHAR(255),
    '''
    os.makedirs(output_directory, exist_ok=True) #exist_ok to avoid errors if the dir already exists
    network_original = pd.read_table(network_table, sep='\t')
    loc_to_at_df = pd.read_table(loc_to_at_file, sep='\t')
    # ------------- genes and homologues -------------
    #gather all genes in my network
    genes1 = set(network_original['Source'].dropna().unique())
    genes2 = set(network_original['Target'].dropna().unique())
    genes_in_network = genes1.union(genes2)
    #Series is used to basicly from a set/list generate a column
    all_genes = sorted(genes_in_network)
    nodes_df = pd.DataFrame({'name': all_genes})
    nodes_df = pd.merge(
                        nodes_df,
                        loc_to_at_df,
                        left_on='name',
                        right_on='Gene_ID',
                        how='left')
    
    # Keep only the 'name' and 'arabidopsis_gene' columns
    nodes_df = nodes_df[['name', 'hit']]
    nodes_df['hit'] = nodes_df['hit'].fillna('--')
    
    # ------------- TF -------------
    tf_df = pd.read_csv(tf_file, header=None, names=['name'])
    nodes_df = pd.merge(
                        nodes_df,
                        tf_df,
                        on='name',
                        how='left',
                        indicator=True) #new column with info about the matches ('both', 'left_only' ...)
    nodes_df['isTF'] = (nodes_df['_merge'] == 'both').astype('Int32')
    #remove the merge column
    nodes_df.drop('_merge', axis=1, inplace=True)

    # ------------- TR -------------
    tr_df = pd.read_csv(tr_file, header=None, names=['name'])
    nodes_df = pd.merge(
                    nodes_df,
                    tr_df,
                    on='name',
                    how='left',
                    indicator=True) #new column with info about the matches ('both', 'left_only' ...)
    
    nodes_df['isTR'] = (nodes_df['_merge'] == 'both').astype('Int32')
    #remove the merge column
    nodes_df.drop('_merge', axis=1, inplace=True)
        
    # ------------- gene_annotation -------------
    gene_annotation_df = pd.read_table(gene_annotation_file, sep='\t')
    nodes_df = pd.merge(
                    nodes_df,
                    gene_annotation_df,
                    left_on='name',
                    right_on='GeneID',
                    how='left')
    
    nodes_df.drop(columns=['GeneID'], inplace=True)
    nodes_df.rename(columns={'Description': 'gene_annotation'}, inplace=True)



    #-----------------
    print()
    print('---- INFO about the NODE DATA ----')
    print()
    len_isTF = (nodes_df['isTF'] == 1).sum()
    print(f'Number TFs genes in the network: {len_isTF}')
    len_isTR = (nodes_df['isTR'] == 1).sum()
    print(f'Number TRs genes in the network: {len_isTR}')


    nodes_df.to_csv(f'{output_directory}/{output_file_name}', index=False, sep=',')



def generate_edges_table(network_table, connecTF_table, cis_table, dap_seq_table, out_dir, output_edges_name):
    """
        I am not using a relational database anymore but this serve as guide to undestand the data anyway
    
        id VARCHAR(255) PRIMARY KEY,
         Source VARCHAR(255),
         Target VARCHAR(255),
         irp_score FLOAT,
         interaction VARCHAR(255),
         width FLOAT,
         connecTF INT,
         edge_type VARCHAR(255),
         gene_name VARCHAR(255),
         cis_elements INT,
         cis_value INT,
         dap_seq INT,
         tf_rank []INT
         directed INT,
    
    """
    pd.set_option('display.max_columns', None)

    network_original = pd.read_table(network_table, sep = '\t')
    network = network_original.iloc[:, :3]
    network['id'] = network['Source'].astype(str) + ' interacts with ' + network['Target'].astype(str)
    network['interaction'] = 'interacts with'
    
# ------------- Connec_TF -------------
    '''  
    connecTF_table = pd.read_table(connecTF_table, sep = '\t')


    # Group by 'gene_id' and 'TARGET' because there were duplicates in this subset
    # - For 'EDGE_TYPE', concatenate the values with a comma
    # - For 'gene_name', keep the first value in each group
    connecTF_table = connecTF_table.groupby(['gene_id', 'TARGET'], as_index=False).agg({
        'EDGE_TYPE': ', '.join,   # Concatenate values in EDGE_TYPE with a comma
        'gene_name': 'first'      # Keep the first value in gene_name
    })

    # Set Pandas option to display all columns
    pd.set_option('display.max_columns', None)

    # Merge for the first case: Source -> gene_id, Target -> TARGET
    merge_1 = pd.merge(network, 
                   connecTF_table,
                   left_on=['Source', 'Target'],
                   right_on=['gene_id', 'TARGET'],
                   how='left')
    
    null_count_col1 = merge_1['gene_id'].isnull().sum()
    #print(merge_1.head(10))

    print(merge_1[merge_1['gene_id'].notnull()])

    # Display the result
    print("Total number of duplicates in merge1:", null_count_col1)

    # Merge for the second case: Source -> TARGET, Target -> gene_id
    merge_2 = pd.merge(network, 
                   connecTF_table,
                   left_on=['Source', 'Target'],
                   right_on=['TARGET', 'gene_id'],
                   how='left')
    
    null_count_col1 = merge_2['gene_id'].isnull().sum()
    #print(merge_2.head(10))

    print('UOUOUOUOUOUOUOOUOU', merge_2[merge_2['gene_id'].notnull()])

    # Display the result
    print("Total number of duplicates in merge1:", null_count_col1)
                  
    print(merge_1.shape)       
    print(merge_2.shape)         
    # Combine both results and remove duplicates (if any)

    common_connecTF_interactions = pd.concat([merge_1, merge_2]).drop_duplicates()
    print(common_connecTF_interactions[common_connecTF_interactions['gene_id'].notnull()].shape)
    print(common_connecTF_interactions.columns)
    print(common_connecTF_interactions.shape)
    # 1. Combine two columns into a set
    #common_connecTF_interactions['connecTF'] = common_connecTF_interactions['gene_id'].isnull().astype('Int32')

    #common_connecTF_interactions = common_connecTF_interactions.groupby(['Source', 'Target'], as_index=False)['connecTF'].sum()
    # print(common_connecTF_interactions.shape)
    # print(common_connecTF_interactions.head())
    # 2. Group by the new column 'combined' and keep the row with the fewest nulls
    def count_nulls(row):
        return row.isnull().sum()
    
    # Find duplicates based on a subset of columns (e.g., 'col1' and 'col2')
    duplicates = common_connecTF_interactions.duplicated(subset=['Source', 'Target'], keep=False)

    # Total number of duplicate rows
    total_duplicates = duplicates.sum()
    print('total duplicates in final')
    print(total_duplicates)

    # Find duplicated rows based on a subset of columns (e.g., 'col1' and 'col2')
    #duplicates = common_connecTF_interactions[common_connecTF_interactions.duplicated(subset=['Source', 'Target'], keep=False)]
    # Identify duplicates based on 'Source' and 'Target'

    # Sort by 'Source' and 'Target'


    # Sort by the count of nulls and keep the row with the fewest nulls - NOT WORKING
    #common_connecTF_interactions = common_connecTF_interactions.loc[common_connecTF_interactions.groupby(['Source', 'Target'])[['irp_score', 'id', 'interaction', 'width', 'EDGE_TYPE', 'gene_name']].apply(count_nulls).idxmin()]
    
    # First, count nulls for each row and create a new column in the dataframe
    common_connecTF_interactions['null_count'] = common_connecTF_interactions[['irp_score', 'id', 'interaction', 'width', 'EDGE_TYPE', 'gene_name']].apply(count_nulls, axis=1)

    
    # Group by 'Source' and 'Target', then get the index of the row with the minimum null count
    idx_min_nulls = common_connecTF_interactions.groupby(['Source', 'Target'])['null_count'].idxmin()

    # Now, filter the rows from the original DataFrame using the idx_min_nulls indices
    common_connecTF_interactions = common_connecTF_interactions.loc[idx_min_nulls]

    print(common_connecTF_interactions.loc[idx_min_nulls].values)
    print('common_connecTF_interactions')
    print(common_connecTF_interactions.shape)
    '''
    connecTF_table = pd.read_table(connecTF_table, sep = '\t')

    merge_1 = pd.merge(network, 
                   connecTF_table,
                   left_on=['Source', 'Target'],
                   right_on=['gene_id', 'TARGET'],
                   how='inner')
    
    merge_2 = pd.merge(network, 
                   connecTF_table,
                   left_on=['Source', 'Target'],
                   right_on=['TARGET', 'gene_id'],
                   how='inner')

    common_connecTF_interactions = pd.concat([merge_1, merge_2]).drop_duplicates()

    # Keep only unique Source-Target pairs
    common_connecTF_interactions = common_connecTF_interactions.drop_duplicates(subset=['Source', 'Target'])
    
    # Create a set of tuples with the common pairs
    common_pairs = set(zip(common_connecTF_interactions['Source'], 
                           common_connecTF_interactions['Target']))
    
    # Check each row in network and create the new column; astype int turns from boolean to 0 and 1
    network['connecTF'] = network[['Source', 'Target']].apply(tuple, axis=1).isin(common_pairs).astype('Int32')
    
    network['EDGE_TYPE'] = '--'
    network['gene_name'] = '--'
    
    # Merge only the rows where connecTF == 1
    merged_data = network.loc[network['connecTF'] == 1, ['Source', 'Target']].merge(
        common_connecTF_interactions, 
        on=['Source', 'Target'], 
        how='left'
    )

    #IF connecTF == 1 im changing the network so that the TF stays in the Source column and the target in the Target. This wont affect the co-exp but will ease knowing the direction of the relation
    network['gene_id'] = '--'  # Default placeholder
    network['TARGET'] = '--'
    # Ensure merged_data length matches network['connecTF'] == 1 length
    if len(merged_data) == len(network.loc[network['connecTF'] == 1]):
        network.loc[network['connecTF'] == 1, ['EDGE_TYPE', 'gene_id', 'TARGET', 'gene_name']] = \
            merged_data[['EDGE_TYPE', 'gene_id', 'TARGET', 'gene_name']].values
        
        # Change source and target
        tf_set = set(connecTF_table['gene_id'])
        # Case 1: both Source and Target are TFs → enforce gene_id → TARGET direction
        mask_both_tf = (network['connecTF'] == 1) & \
                    (network['Source'].isin(tf_set)) & (network['Target'].isin(tf_set))

        network.loc[mask_both_tf, 'Source'] = network.loc[mask_both_tf, 'gene_id']
        network.loc[mask_both_tf, 'Target'] = network.loc[mask_both_tf, 'TARGET']

        # Case 2: only Target is TF → flip direction
        mask_target_tf = (network['connecTF'] == 1) & \
                        (network['Target'].isin(tf_set)) & (~network['Source'].isin(tf_set))

        network.loc[mask_target_tf, ['Source', 'Target']] = network.loc[mask_target_tf, ['Target', 'Source']].values

    else:
        print("Warning: Mismatch in row counts, possible issue with merge.")    
  
    network = network.drop(columns=['gene_id', 'TARGET'])

    
    len_connecTF = (network['connecTF'] == 1).sum()

    # ------------- cis_elements -------------
    
    network['cis_elements'] = 0
    network['cis_value'] = 0

    cis_table = pd.read_table(cis_table, sep='\t', header=None)
    cis_table.columns = ['gene_id', 'TARGET', 'cis_bool']

    merge_1 = pd.merge(network, 
                   cis_table.drop_duplicates(),
                   left_on=['Source', 'Target'],
                   right_on=['gene_id', 'TARGET'])

    # Merge for the second case: Source -> TARGET, Target -> gene_id
    merge_2 = pd.merge(network, 
                   cis_table.drop_duplicates(),
                   left_on=['Source', 'Target'],
                   right_on=['TARGET', 'gene_id'])
    
    common_cis_interactions = pd.concat([merge_1, merge_2]).drop_duplicates()    
   
    common_cis_interactions['cis_value'] = common_cis_interactions['cis_bool'].apply(lambda x: float(x.split('_')[1]) if x != 'NO' else 0)
    common_cis_interactions['cis_bool'] = common_cis_interactions['cis_bool'].apply(lambda x: 1 if x != 'NO' else 0)

    # Update cis_elements based on matching pairs in common_cis_interactions
    merged_cis = network[['Source', 'Target']].merge(common_cis_interactions[['Source', 'Target', 'cis_bool', 'cis_value']], 
                                                    on=['Source', 'Target'], how='left')

    network['cis_elements'] = merged_cis['cis_bool'].fillna(0)  # Default to 0 if no match is found
    #turn the data to int to get 1 ontead of 1.0
    network['cis_elements'] = network['cis_elements'].astype('Int32')

    #Flip direction if TF is in Target and there is no ConnecTF evidence

    network['gene_id'] = '--'
    network['TARGET'] = '--'

    # Prepare values to fill in
    if len(common_cis_interactions) >= len(network.loc[network['cis_elements'] == 1]):
        merged_vals = network.loc[network['cis_elements'] == 1, ['Source', 'Target']].merge(
            common_cis_interactions[['Source', 'Target', 'gene_id', 'TARGET']],
            on=['Source', 'Target'],
            how='left'
        )
        network.loc[network['cis_elements'] == 1, ['gene_id', 'TARGET']] = merged_vals[['gene_id', 'TARGET']].values

        # Identify rows to flip where no ConnecTF evidence exists
        mask = (
            ((network['Source'] == network['TARGET']) | (network['Target'] == network['gene_id'])) &
            (network['connecTF'] == 0)
        )

        # Flip Source/Target only for rows matching the mask
        network.loc[mask, ['Source', 'TARGET']] = network.loc[mask, ['TARGET', 'Source']].values
        network.loc[mask, ['Target', 'gene_id']] = network.loc[mask, ['gene_id', 'Target']].values
    else:
        print("Warning: Mismatch in cis-element merge counts.")

    # Drop temporary columns
    network = network.drop(columns=['gene_id', 'TARGET'])
    
    # ------------- DAP Seq data -------------
    dap_seq_df = pd.read_table(dap_seq_table, sep='\t')
    network['dap_seq'] = 0

    merged_net_dap = pd.merge(network, 
                   dap_seq_df,
                   left_on=['Source', 'Target'],
                   right_on=['gene_id', 'TARGET'])
    merged_net_dap['dap_bool'] = 1
    # Update dap_seq based on matching pairs in common_dap_interactions
    merged_dap = network[['Source', 'Target']].merge(
                    merged_net_dap[['Source', 'Target', 'dap_bool']], 
                    on=['Source', 'Target'], 
                    how='left')

    network['dap_seq'] = merged_dap['dap_bool'].fillna(0).astype('Int32')

    # List of dap-seq TFs
    tfs = ['LOC111997151', 'LOC112008346', 'LOC112030452']

    # Identify rows with DAP-seq evidence
    dap_supported = network['dap_seq'] == 1

    # Identify rows where the TF is in the Target
    tf_in_target = dap_supported & network['Target'].isin(tfs) & ~network['Source'].isin(tfs)

    # Swap Source and Target for those rows
    network.loc[tf_in_target, ['Source', 'Target']] = network.loc[tf_in_target, ['Target', 'Source']].values



    ## ------------- TF_Rank data -------------

    network['tf_rank'] = network.apply(determine_tf_rank, axis=1)
    
    ## ------------- Directed -------------
    
    network['directed'] = network.apply(determine_direction, axis=1)
    print(network[network['dap_seq'] == 1][['Source', 'Target', 'tf_rank']])

    print()
    print('---- INFO about the EDGE DATA ----')
    print()
    filtered_rows = network[network['tf_rank'].apply(lambda x: 3 in x)]
    print(filtered_rows)
    #print('Edge Table:')
    #print(network.head())

    len_connecTF = (network['connecTF'] == 1).sum()
    print(f'Number of predicted connecTF TFs interactions in the network: {len_connecTF}')
    # If matching rows exist, print them
    
    print()
    len_cis = (network['cis_elements'] == 1).sum()
    print(f'Number of predicted TFs interactions through cis-elements analysis in the network: {len_cis}')
    print()

    len_dap = (( network['dap_seq'] == 1).sum())
    print(f'Number of  TFs interactions achieved by Dap Seq technique in the network: {len_dap}')
    
    len_dir = (network['directed'] == 1).sum()
    print(f'Number of directed edges in the network: {len_dir}')
    #write table to file
    network.to_csv(f'{out_dir}/{output_edges_name}', index=False, sep=',')
    print()
    print(f"edges.csv saved on {out_dir}")

    

import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Script that builds the nodes and edges tables for the network.",
        epilog="""
                Examples:

                nodes table:
                    python generate_edges_nodes_tables.py --function generate_nodes_table 
                        --network data/seidr_output/network.csv \\
                        --tf_file data/connecTF_files/tf_file.txt \\
                        --tr_file data/connecTF_files/tr_file.txt \\
                        --gene_annotation_file data/genome/genome_annotation/gene_functional_annotation.csv \\
                        --loc_to_at_file data/connecTF_files/LOC_to_Arabidopsis_match.txt \\
                        --out_dir output \\
                        --output_nodes_name

                edges table:
                    python generate_edges_nodes_tables.py --function generate_edges_table 
                        --network data/seidr_output/network.csv \\
                        --connecTF_table data/connecTF_files/connecTF_results.txt \\
                        --cis_table data/connecTF_files/promoter_motif_enrichment.txt \\
                        --dap_seq_table data/connecTF_files/dap_seq_results.txt \\
                        --out_dir output \\
                        --output_edges_name

                Use --help with any argument for more info.
                """
    )
    parser.add_argument('--function', required=True, choices=[
        'generate_nodes_table',
        'generate_edges_table',
    ], help='Function to run.')
    
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    # Global options
    parser.add_argument('--network', type=str, default=f'{os.path.join(base_dir, 'seidr_output', 'network.csv')}',help='Path to network file.')
    parser.add_argument('--out_dir', type=str, default=f'{os.path.join(base_dir, 'connecTF','connecTF_Targeted', 'latest', 'LOC_to_Arabidopsis_match.txt')}', help='Path to file that has a columns with locs and another with the corresponding arabidopsis homologues.')

    # generate_nodes_table options
    parser.add_argument('--tf_file', type=str, default=f'{os.path.join(base_dir, 'connecTF','connecTF_Targeted', 'latest', 'tf_file.txt')}', help='Path to the file with the nodes that are tfs.')
    parser.add_argument('--tr_file', type=str, default=f'{os.path.join(base_dir, 'connecTF','connecTF_Targeted', 'latest', 'tr_file.txt')}', help='Path to the file with the nodes that are trs.')
    parser.add_argument('--gene_annotation_file', type=str, default=f'{os.path.join(base_dir, 'genome', 'genome_annotation', 'gene_functional_annotation.csv')}', help='Path to script with the nodes annotations previously extracted from the quercus suber genome GTF file.')
    parser.add_argument('--loc_to_at_file', type=str, default=f'{os.path.join(base_dir, 'connecTF','connecTF_Targeted', 'latest', 'LOC_to_Arabidopsis_match.txt')}', help='Path to file that has a columns with locs and another with the corresponding arabidopsis homologues.')
    parser.add_argument('--output_nodes_name', type=str, default='nodes.csv', help='Name of the output nodes file.')

    # generate_edges_table options
    parser.add_argument('--connecTF_table', type=str, default=f'{os.path.join(base_dir, 'connecTF', 'connecTF_Targeted','latest', 'output_targeted', 'QS_output_final.csv')}', help='Path to connecTF output table.')
    parser.add_argument('--cis_table', type=str, default=f'{os.path.join(base_dir, 'cis_elements', 'fimo_out', 'ALL_TF_Matches_Output.txt')}', help='Path to cis-elements analysis output table.')
    parser.add_argument('--dap_seq_table', type=str, default=f'{os.path.join(base_dir, 'dap_seq', 'dap_seq.txt')}',help='Path to DAPseq output file.')
    parser.add_argument('--output_edges_name', type=str, default='edges.csv', help='Name of the output edges file.')


    args = parser.parse_args()


    if args.function == "generate_nodes_table":
        generate_nodes_table(
            args.network,
            args.tf_file,
            args.tr_file,
            args.gene_annotation_file,
            args.loc_to_at_file,
            args.out_dir,
            args.output_nodes_name
        )
    elif args.function == "generate_edges_table":
        generate_edges_table(
            args.network,
            args.connecTF_table,
            args.cis_table,
            args.dap_seq_table,
            args.out_dir,
            args.output_edges_name
        )
    else:
        parser.print_help()

        
if __name__ == "__main__":
    main()
