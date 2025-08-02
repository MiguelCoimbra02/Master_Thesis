#!/bin/bash
#This script converts SRA files into FASTQ format using the fasterq-dump tool from NCBI's SRA Toolkit. 
#It processes each SRR from a given BioProject, converting them into FASTQ files.

# ./fasterq_dump.sh Bioproject_accession
#------------------------------------------------------------------------------
BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BIN_DIR="$(cd "$(dirname "$0")/../../bin" && pwd)"


# Path to the program fasterq-dump
SRA_TOOLKIT_DIR="$BASE_DIR/../bin/sratoolkit.3.1.1-ubuntu64/bin"
export PATH=$PATH:"$BIN_DIR/sratoolkit.3.1.1-ubuntu64/bin"

# Base directory for BioProjects

# Check if the BioProject name is provided as an argument
if [ -z "$1" ]; then
    # If not provided, prompt the user
    read -p "Enter the BioProject name: " BIOPROJECT
else
    # Use the argument as the BioProject name
    BIOPROJECT="$1"
fi

BIOPROJECT_DIR="$BASE_DIR/$BIOPROJECT"

# Check if the BioProject directory exists
if [ ! -d "$BIOPROJECT_DIR" ]; then
    echo "Error: BioProject directory '$BIOPROJECT' not found. Run the command ./prefetch_sra.sh <$BIOPROJECT> before."
    exit 1
fi

SRR_LIST="$BASE_DIR/bioProjects_info/$BIOPROJECT.txt"

if [ ! -f "$SRR_LIST" ]; then
    echo "Error: SRR list file '$SRR_LIST' not found."
    exit 1
fi

# Ensure the output directory exists
mkdir -p "$BIOPROJECT_DIR/fastq"

# Loop through each SRR in raw_samples.txt
while IFS= read -r SRR; do
    echo "Processing SRR: $SRR for BioProject: $BIOPROJECT"
    
    # Run fasterq-dump for the current SRR
    # Use absolute path for clarity
    if ! fasterq-dump "$BASE_DIR/$BIOPROJECT/prefetch/$SRR" -O "$BIOPROJECT_DIR/fastq"; then
        echo "Error: fasterq-dump failed for $SRR"
        continue
    fi
done < "$SRR_LIST"

echo "Fasterq-dump complete for BioProject: $BIOPROJECT."
