#!/bin/bash
#This script processes a featureCounts summary file to extract the number of assigned and unassigned reads, adjusts the values for paired-end (PE) reads
#and logs the results into a CSV file for a specified BioProject.

# ./log_fc.sh /path/to/summary.txt format Bioproject
#------------------------------------------------------------------------------

RAW_READS_SUMMARY="$1"
FORMAT="$2"
BIOPROJECT="$3"

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

if [ "$FORMAT" == "pe" ]; then
    # Write to CSV file for paired-end reads
    printf "Assigned Reads (pair) \n" > $BASE_DIR/$BIOPROJECT/logs/a.csv
    printf "Unassigned Reads (No Features) \n" > $BASE_DIR/$BIOPROJECT/logs/u.csv

else
    # Write to CSV file for single-end reads
    printf "Assigned Reads \n" > $BASE_DIR/$BIOPROJECT/logs/a.csv
    printf "Unassigned Reads (No Features) \n" > $BASE_DIR/$BIOPROJECT/logs/u.csv

fi

n_assigned_r=$(grep "Assigned" $RAW_READS_SUMMARY | grep -o '[[:digit:]]*')
n_unassigned_r=$(grep "Unassigned_NoFeatures" $RAW_READS_SUMMARY | grep -o '[[:digit:]]*')

echo "$n_assigned_r" >> $BASE_DIR/$BIOPROJECT/logs/a.csv
echo "$n_unassigned_r" >> $BASE_DIR/$BIOPROJECT/logs/u.csv


paste $BASE_DIR/$BIOPROJECT/logs/a.csv $BASE_DIR/$BIOPROJECT/logs/u.csv > $BASE_DIR/$BIOPROJECT/logs/featureCounts_log.csv

rm $BASE_DIR/$BIOPROJECT/logs/a.csv
rm $BASE_DIR/$BIOPROJECT/logs/u.csv