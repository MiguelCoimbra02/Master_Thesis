#!/bin/bash
#This script consolidates logs from various RNA-seq processing steps into a single CSV file for a given BioProject.

# ./log_BioP_stats.sh Bioproject
#------------------------------------------------------------------------
BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

bioP="$1"
bioP_dir="$BASE_DIR/$bioP/logs"
paste $bioP_dir/STAR_log.csv $bioP_dir/featureCounts_log.csv > $bioP_dir/${bioP}_Statistics.tsv
#rm $bioP_dir/trimmomatic_log.csv
#rm $bioP_dir/STAR_log.csv
#rm $bioP_dir/featureCounts_log.csv