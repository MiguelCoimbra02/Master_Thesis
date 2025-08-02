#!/bin/bash
#Parses a STAR log file to extract alignment metrics such as uniquely mapped reads and unmapped reads.
#Then logs these metrics into a CSV file for a given BioProject.

# ./log_star.sh path/to/sample.Log.final.out BioProject format(se;pe)
#------------------------------------------------------------------------------
# Extract sample name
SRR="$(basename "$1" .Log.final.out)"
echo "Processing sample: $SRR"
star_log_file="$1"
BIOPROJECT="$2"

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
FORMAT="$3"



input_reads_n=$(grep "Number of input reads" $star_log_file | grep -o '[[:digit:]]*')
uniq_map_reads_n=$(grep "Uniquely mapped reads number" $star_log_file | grep -o '[[:digit:]]*')
uniq_map_reads_p=$(grep "Uniquely mapped reads %" $star_log_file | grep -o '[[:digit:]]\+.[[:digit:]]*')
multi_mapped_reads_n=$(grep "Number of reads mapped to multiple loci" $star_log_file | grep -o '[[:digit:]]*')
multi_mapped_reads_p=$(grep "% of reads mapped to multiple loci" $star_log_file | grep -o '[[:digit:]]\+.[[:digit:]]*')
unmapped_reads_short=$(grep "Number of reads unmapped: too short" $star_log_file | grep -o '[[:digit:]]*')
unmapped_reads_other=$(grep "Number of reads unmapped: other" $star_log_file | grep -o '[[:digit:]]*')
mismatch_rate=$(grep "Mismatch rate per base, %" $star_log_file | grep -o '[[:digit:]]\+.[[:digit:]]*')

# Append data to the CSV file
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$SRR" "$FORMAT" "$input_reads_n" "$uniq_map_reads_n" "$uniq_map_reads_p" "$multi_mapped_reads_n" "$multi_mapped_reads_p" "$unmapped_reads_short" "$unmapped_reads_other" "$mismatch_rate" >> $BASE_DIR/$BIOPROJECT/logs/STAR_log.csv
