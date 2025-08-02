#!/bin/bash
# This script fetches SRA files for a given BioProject using the prefetch tool from NCBI's SRA Toolkit.
# It runs the prefetch command in parallel for each SRR in the specified BioProject and logs the process.

# ./prefetch.sh Bioproject_accession

#------------------------------------------------------------------------------ 

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BIN_DIR="$(cd "$(dirname "$0")/../../bin" && pwd)"

# Path to the program prefetch
export PATH="$BIN_DIR/sratoolkit.3.1.1-ubuntu64/bin:$PATH"

# Check if the BioProject name is provided as an argument
if [ -z "$1" ]; then
    # If not provided, prompt the user
    read -p "Enter the BioProject name: " BIOPROJECT
else
    # Use the argument as the BioProject name
    BIOPROJECT="$1"
fi

# Define the file that contains SRR numbers for the given BioProject
SRA_FILE="$BASE_DIR/bioProjects_info/${BIOPROJECT}.txt"

# Ensure that the SRA file exists
if [ ! -f "$SRA_FILE" ]; then
    echo "SRA file for BioProject '$BIOPROJECT' not found!"
    exit 1
fi

# Create the main directory with the BioProject name
MAIN_DIR="$BASE_DIR/$BIOPROJECT"
mkdir -p "$MAIN_DIR/logs"
PREFETCH_DIR="$MAIN_DIR/prefetch"

# Ensure the directories exist
mkdir -p "$PREFETCH_DIR"

# Define the log file path
LOG_FILE="$MAIN_DIR/logs/prefetch_log.txt"

# Count the number of SRR entries in the SRA file and set the number of parallel jobs
PARALLEL_JOBS=$(wc -l < "$SRA_FILE")
if [ "$PARALLEL_JOBS" -gt 14 ]; then
    PARALLEL_JOBS=14
fi

echo "Number of parallel jobs set to: $PARALLEL_JOBS"

# Function to prefetch only if the directory is empty
prefetch_if_needed() {
    SRR="$1"
    SRR_DIR="$PREFETCH_DIR/$SRR"
    
    if [ -d "$SRR_DIR" ] && [ "$(ls -A "$SRR_DIR" 2>/dev/null)" ]; then
        echo "$SRR directory is not empty. Skipping prefetch."
    else
        prefetch --max-size 100G --output-directory "$PREFETCH_DIR" "$SRR" && echo "$SRR successfully prefetched." || echo "Failed to prefetch $SRR."
    fi
}

export -f prefetch_if_needed
export PREFETCH_DIR

# Run prefetch in parallel, skipping if the SRR directory is not empty
cat "$SRA_FILE" | parallel -j "$PARALLEL_JOBS" prefetch_if_needed {} >> "$LOG_FILE" 2>&1

# Notify the user when done
echo "Prefetching completed. Check the log file at '$LOG_FILE' for details."
