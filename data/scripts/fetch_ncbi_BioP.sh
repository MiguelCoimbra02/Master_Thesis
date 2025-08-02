#!/bin/bash

#This script fetches the BioProject IDs and related SRA run data from NCBI based on a search query. 
#It saves a list of BioProject names to a text file and creates a summary file with BioProject titles and SRA run information (e.g., paired-end or single-end).
#All the following Scripts will handle both single-end and paired-end reads

# ./fetch_ncbi_BioP.sh or ./fetch_ncbi_BioP.sh 'Custom NCBI Query'
#------------------------------------------------------------------------------

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"


# Define the default search query
DEFAULT_QUERY='("Quercus suber"[Organism] AND "biomol rna"[Properties] AND ("platform bgiseq"[Properties] OR "platform illumina"[All Fields])'

# Check if a query argument is passed to the script, else use the default
if [ $# -gt 0 ]; then
    QUERY="$1"
    echo "Using custom query: $QUERY"
else
    QUERY="$DEFAULT_QUERY"
    echo "Using default query."
fi
# Define the output file path
mkdir -p $BASE_DIR/bioProjects_info
OUTPUT_FILE_list="$BASE_DIR/bioProjects_info/COak_BioP_list.txt"
OUTPUT_FILE_summary="$BASE_DIR/bioProjects_info/COak_BioP_summary_unordered.txt"


# Check if the output file already exists, if not, fetch and save the bioProjects names
if [ ! -f "$OUTPUT_FILE_list" ]; then
    # Run esearch and efetch to get bioProjects run IDs and save them to a text file
    esearch -db sra -query "$QUERY" | efetch -format runinfo | tail -n +2 | cut -d "," -f22 | sort | uniq > "$OUTPUT_FILE_list"
    
    # Notify the user that the SRA run IDs have been saved
    echo "SRA BioProject names have been saved to $OUTPUT_FILE_list"
else
    # Notify the user that the file already exists and the write is skipped
    echo "File already exists, skipping write."
fi

# Create or overwrite the summary file
> "$OUTPUT_FILE_summary"

# Read BioProjects from the file and iterate with a `for` loop
for BIOPROJECT in $(cat "$OUTPUT_FILE_list"); do
    if [ -n "$BIOPROJECT" ]; then

        # Fetch SRR data for the current BioProject
        BIOPROJECT_TITLE=$(esearch -db bioproject -query "${BIOPROJECT}[BioProject]" | efetch -format docsum | xtract -pattern DocumentSummary -element Title)
        SRR_DATA=$(esearch -db sra -query "${BIOPROJECT}[BioProject]" | efetch -format runinfo | tail -n +2 | awk -F ',' '{print $1, $16}')

        # Write the BioProject name, title, and date to the summary file
        echo "BioProject: $BIOPROJECT" >> "$OUTPUT_FILE_summary"
        echo "Title: $BIOPROJECT_TITLE" >> "$OUTPUT_FILE_summary"

        # Create a new file named after the BioProject (e.g., PRJNA347903.txt)
        BIOPROJECT_FILE="$BASE_DIR/bioProjects_info/${BIOPROJECT}_unordered.txt"
        

        # Check if any SRRs were found for the current BioProject
        if [ -n "$SRR_DATA" ]; then
            # Process each line of SRR data
            echo "$SRR_DATA" | while read -r SRR STATUS; do
                # Write each SRR to the BioProject-specific file
                echo "$SRR" >> "$BIOPROJECT_FILE"

                # Write to the summary file with (PE) or (SE)
                if [ "$STATUS" == "PAIRED" ]; then
                    echo "--> $SRR (PE)" >> "$OUTPUT_FILE_summary"
                else
                    echo "--> $SRR (SE)" >> "$OUTPUT_FILE_summary"
                fi
            done

            echo "" >> "$OUTPUT_FILE_summary" # Add a blank line for separation
            sort -n "$OUTPUT_FILE_summary" > "$BASE_DIR/bioProjects_info/COak_BioP_summary.txt"
            rm "$OUTPUT_FILE_summary"
        else
            # Debug: No SRRs found for the current BioProject
            echo "No SRRs found for BioProject: $BIOPROJECT"
        fi
    else
        # Debug: Empty line in the BioProject list
        echo "Skipping empty line in BioProject list."
    fi
done

# Notify the user that the summary has been created
echo "BioProject(s) summary has been saved to $BASE_DIR/bioProjects_info"