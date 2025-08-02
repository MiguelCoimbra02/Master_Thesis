#!/bin/bash

# Runs a sequence of scripts for a given BioProject: prefetch, fasterq_dump, preprocessing, STAR alignment, and cleaning.
# Usage: ./bioProject_pipeline.sh BioProject [format]
#------------------------------------------------------------------------------

if [ -z "$1" ]; then
    read -p "Enter the BioProject name: " BIOPROJECT
else
    BIOPROJECT="$1"
fi

if [ -z "$2" ]; then
    # If not provided, prompt the user
    read -p "Enter the format of the reads (se or pe) for BioProject $BIOPROJECT: " FORMAT
else
    # Use the argument as the format
    FORMAT="$2"
fi

SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"

SCRIPTS=(
    "$SCRIPTS_DIR/prefetch.sh $BIOPROJECT"
    "$SCRIPTS_DIR/fasterq_dump.sh $BIOPROJECT"
    "$SCRIPTS_DIR/preprocessing.sh $BIOPROJECT $FORMAT"
    "$SCRIPTS_DIR/star_fc.sh $BIOPROJECT $FORMAT"
    "$SCRIPTS_DIR/clean_BioP.sh $BIOPROJECT"
)

for script_command in "${SCRIPTS[@]}"; do
    script_name=$(echo "$script_command" | awk '{print $1}')
    if [ -x "$script_name" ]; then
        echo "Running $script_command..."
        # Run and wait for the script to complete
        eval "$script_command" || { echo "Error running $script_name"; exit 1; }
    else
        echo "Script $script_name is not executable. Please check permissions."
        exit 1
    fi
done
echo "All scripts completed successfully."