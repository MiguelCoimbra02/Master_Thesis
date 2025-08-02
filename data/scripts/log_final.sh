#!/bin/bash
#This script generates a Star Log. Created to be run alone, not part of other script. 
#It greps information from the .Log.final.out (STAR)

# ./log_final.sh Bioproject_accession format (se; pe)
#------------------------------------------------------------------------------
# The input STAR log file
BIOPROJECT=$1  # Takes the first argument passed to the script as the log file
FORMAT=$2
BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SRR_LIST="$BASE_DIR/bioProjects_info/$BIOPROJECT.txt"

STAR_LOG_DIR="$BASE_DIR/$BIOPROJECT/aligned_reads"
FC_LOG_DIR="$BASE_DIR/$BIOPROJECT/logs"
OUT_CSV="$FC_LOG_DIR"



if [ "$FORMAT" == "pe" ]; then
    # Write to CSV file for paired-end reads
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "Sample Name" "Format" "Total Reads (pair)" "Uniquely Mapped Reads" "Uniquely Mapped %" \
    "Multi-Mapped Reads" "Multi-Mapped %" "Unmapped Reads (Too Short)" "Unmapped Reads (Other)" \
    "Mismatch Rate (%)" > $OUT_CSV/final_log_$BIOPROJECT.csv

    printf "Assigned Reads (pair) \n" > $BASE_DIR/$BIOPROJECT/logs/a.csv
    printf "Unassigned Reads (No Features) \n" > $BASE_DIR/$BIOPROJECT/logs/u.csv

else
    # Write to CSV file for single-end reads
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "Sample Name" "Format" "Total Reads" "Uniquely Mapped Reads" "Uniquely Mapped %" \
    "Multi-Mapped Reads" "Multi-Mapped %" "Unmapped Reads (Too Short)" "Unmapped Reads (Other)" \
    "Mismatch Rate (%)" > $OUT_CSV/final_log_$BIOPROJECT.csv

    printf "Assigned Reads \n" > $BASE_DIR/$BIOPROJECT/logs/a.csv
    printf "Unassigned Reads (No Features) \n" > $BASE_DIR/$BIOPROJECT/logs/u.csv

fi


# Read SRR list and process each entry
while IFS= read -r SRR; do
    star_log_file="$STAR_LOG_DIR/"$SRR".Log.final.out"
    
    # Retrieve information from the STAR log file using grep
    input_reads_n=$(grep "Number of input reads" $star_log_file | grep -o '[[:digit:]]*')
    uniq_map_reads_n=$(grep "Uniquely mapped reads number" $star_log_file | grep -o '[[:digit:]]*')
    uniq_map_reads_p=$(grep "Uniquely mapped reads %" $star_log_file | grep -o '[[:digit:]]\+.[[:digit:]]*')
    multi_mapped_reads_n=$(grep "Number of reads mapped to multiple loci" $star_log_file | grep -o '[[:digit:]]*')
    multi_mapped_reads_p=$(grep "% of reads mapped to multiple loci" $star_log_file | grep -o '[[:digit:]]\+.[[:digit:]]*')
    unmapped_reads_short=$(grep "Number of reads unmapped: too short" $star_log_file | grep -o '[[:digit:]]*')
    unmapped_reads_other=$(grep "Number of reads unmapped: other" $star_log_file | grep -o '[[:digit:]]*')
    mismatch_rate=$(grep "Mismatch rate per base, %" $star_log_file | grep -o '[[:digit:]]\+.[[:digit:]]*')

    # Append data to the CSV file
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$SRR" "$FORMAT" "$input_reads_n" "$uniq_map_reads_n" "$uniq_map_reads_p" "$multi_mapped_reads_n" "$multi_mapped_reads_p" "$unmapped_reads_short" "$unmapped_reads_other" "$mismatch_rate" >> $OUT_CSV/final_log_$BIOPROJECT.csv

done < "$SRR_LIST"
RAW_READS_SUMMARY="$FC_LOG_DIR/"$BIOPROJECT"_raw_counts.txt.summary"
n_assigned_r=$(grep "Assigned" $RAW_READS_SUMMARY | grep -o '[[:digit:]]*')
n_unassigned_r=$(grep "Unassigned_NoFeatures" $RAW_READS_SUMMARY | grep -o '[[:digit:]]*')

echo "$n_assigned_r" >> $BASE_DIR/$BIOPROJECT/logs/a.csv
echo "$n_unassigned_r" >> $BASE_DIR/$BIOPROJECT/logs/u.csv
paste $BASE_DIR/$BIOPROJECT/logs/a.csv $BASE_DIR/$BIOPROJECT/logs/u.csv > $BASE_DIR/$BIOPROJECT/logs/featureCounts_log.csv
paste $OUT_CSV/final_log_$BIOPROJECT.csv $BASE_DIR/$BIOPROJECT/logs/featureCounts_log.csv > $OUT_CSV/final_log_complete_$BIOPROJECT.tsv


echo "Results written to $BIOPROJECT.csv placed in $OUT_CSV."

rm $BASE_DIR/$BIOPROJECT/logs/a.csv
rm $BASE_DIR/$BIOPROJECT/logs/u.csv
rm $OUT_CSV/final_log_$BIOPROJECT.csv
rm $BASE_DIR/$BIOPROJECT/logs/featureCounts_log.csv
