#!/bin/bash
# Fetches a motif matrix for each TF in AT format on PlantTFDB
# Run: get_tf_motif_matrix.sh $1 - TF file with one per line in AT format $2 - a file where i have 2 columns: loc | at match
#-----
#./get_tf_motif_matrix.sh ../connecTF/TF_AT.txt ../connecTF/LOC_to_Arabidopsis.txt

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$BASE_DIR/cis_elements/TF_motifs"

# Check if the user provided an input file
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 TF_AT_FILE loc_at_match.txt"
    echo "Please provide the paths to both the TF list file and the LOC to Arabidopsis match file."
    exit 1
fi

TF_AT_FILE=$1
LOC_AT_MATCH_FILE=$2

# Read the LOC to Arabidopsis gene match file into an associative array
declare -A loc_to_at_map

# Populate the loc_to_at_map to handle cases where multiple LOCs map to the same Arabidopsis gene
while IFS=$'\t' read -r loc at_gene; do
    loc_to_at_map["$at_gene"]+="$loc "
done < "$LOC_AT_MATCH_FILE"

# Loop over each TF in the TF_AT_FILE
while read -r TF; do
    echo "Requesting TF motifs for $TF from PlantTFDB"
    # Download the TF motif
    #wget -P "$BASE_DIR/cis_elements/TF_motifs" "http://planttfdb.gao-lab.org/motif/Ath/$TF.meme"

    # If the TF motif was successfully downloaded, proceed with renaming based on LOC
    if [ -f "$BASE_DIR/cis_elements/TF_motifs/$TF.meme" ]; then
        # Check if the TF matches any Arabidopsis gene
        if [ -n "${loc_to_at_map[$TF]}" ]; then
            # Loop through all LOCs matching this Arabidopsis gene and create a separate file for each LOC
            for loc in ${loc_to_at_map[$TF]}; do
                echo "Match found: $TF -> $loc"
                cp "$BASE_DIR/cis_elements/TF_motifs/$TF.meme" "$BASE_DIR/cis_elements/TF_motifs/$loc.meme"
                #rm "$BASE_DIR/cis_elements/TF_motifs/$TF.meme" 
                echo "Copied $TF.meme to $loc.meme"
            done
        fi
    fi
done < "$TF_AT_FILE"
