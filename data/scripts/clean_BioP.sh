#!/bin/bash

#!/bin/bash
# Removes unnecessary files and directories from a given BioProject directory.
# Usage: ./clean_BioP.sh BioProject
#------------------------------------------------------------------------------

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
 

# Check if the BioProject name is provided as an argument
if [ -z "$1" ]; then
    # If not provided, prompt the user
    read -p "Enter the BioProject name: " BIOPROJECT
else
    # Use the argument as the BioProject name
    BIOPROJECT="$1"
fi

# Check if the BioProject directory exists
BIOPROJECT_DIR="$BASE_DIR/$BIOPROJECT"

if [ ! -d "$BIOPROJECT_DIR" ]; then
    echo "Error: BioProject directory '$BIOPROJECT' does not exist inside '$BASE_DIR'. Please check the directory name."
    exit 1
fi

rm -rf $BIOPROJECT_DIR/fastq
rm -rf $BIOPROJECT_DIR/prefetch
rm -rf $BIOPROJECT_DIR/preprocessing

echo "All the unnecessary files and directories of Bioproject $BIOPROJECT were successfully removed."


