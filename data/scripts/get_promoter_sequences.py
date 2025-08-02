#!/usr/bin/python3
import os
import pandas as pd
# Get the parent directory of the current script
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
out_dir = os.path.join(base_dir, 'cis_elements', 'promoter_sequences.csv')
promoter_sequences_file = os.path.join(base_dir, 'cis_elements', 'promoters', 'promotor_2000_sequences.bed')
promoter_cordinates_file = os.path.join(base_dir, 'cis_elements', 'promoters', 'promotor_2000_pb.bed')
#network_file = os.path.join(base_dir, 'seidr_output', 'network.txt')

def generate_promoter_sequences(promoter_sequences_file = promoter_sequences_file, promoter_cordinates_file = promoter_cordinates_file, output_dir = out_dir):

    """
    Generates a DataFrame containing promoter sequences linked to their corresponding genes.

    The function processes promoter sequence files, extracts scaffold and sequence information, 
    and merges it with gene annotations. The output DataFrame is saved as a CSV file.

    Arguments:
    promoter_sequences_file (str) : Path to the file containing promoter sequences in .bed format.
    promoter_cordinates_file (str) : Path to the file containing promoter coordinates and gene annotations in .bed format.
    output_dir (str) : Directory where the resulting DataFrame will be saved as a CSV file.

    Requires:
    - The input files (`promoter_sequences_file`and `promoter_cordinates_file`) must exist at the specified locations.
    - The promoter sequence file (`promoter_sequences_file`) must be in a BED format with sequences in alternating rows.
    - The promoter coordinates file (`promoter_cordinates_file`) must include columns: `Scaffold`, `Start`, `End`, `Gene`, and `Strand`.
    
    Ensures:
    - The output DataFrame (`query_promoter_sequences`) is saved as a CSV file at the path defined by `output_dir`.
    - The DataFrame contains columns: `Sequence`, `Scaffold`, and `Gene` after merging and processing the input files.
    - The file is saved as `Promoter_sequences_df.csv` within the directory `output_dir`.
    - All operations will be performed in memory, and the function does not return anything.

    Example:
    >>>generate_promoter_sequences(promoter_sequences_file="promoter_sequences.bed", 
                                     promoter_cordinates_file="promoter_coordinates.bed",
                                     network_file="network.txt",
                                     output_dir="/path/to/output")

    The resulting CSV file will contain a merged dataset of promoter sequences and gene annotations.

    """
    #network = pd.read_table(network_file)
    promotor_sequences = pd.read_table(promoter_sequences_file, sep="\t", names=["Scaffolds"])
    promoter_cordinates = pd.read_table(promoter_cordinates_file, sep="\t", names=["Scaffold","Start","End","Gene","Strand"])

    #Get scafold rows
    scaffolds = promotor_sequences.loc[::2].copy()  # Make an explicit copy
    scaffolds.reset_index(drop=True, inplace=True)
    #removes the '>'
    scaffolds["New_Scaffolds"] = scaffolds["Scaffolds"].str[1:]
    #deletes the original collumn with >
    scaffolds = scaffolds.drop("Scaffolds", axis=1)

    sequences = promotor_sequences.loc[1::2]
    sequences.reset_index(drop = True, inplace = True)

    new_promoter_sequences = pd.merge(scaffolds, sequences, left_index = True, right_index = True)
    new_promoter_sequences.columns = ["Scaffold","Sequence"]

    nice_promoter_sequences = pd.concat([new_promoter_sequences, new_promoter_sequences["Scaffold"].str.split(':', expand = True)], axis = 1)
    nice_promoter_sequences = nice_promoter_sequences.drop(nice_promoter_sequences.columns[[0,3]], axis = 1)
    nice_promoter_sequences.columns = ["Sequence", "Scaffold"]

    scaffold_gene_link = promoter_cordinates[["Scaffold", "Gene"]]
    genes = scaffold_gene_link["Gene"]
    query_promoter_sequences = pd.merge(nice_promoter_sequences, genes, left_index = True, right_index = True)
    query_promoter_sequences.to_csv(output_dir, index=False)    

#Run the function to get the tables
generate_promoter_sequences()