#!/bin/bash
#Analyzes a Trimmomatic log file to calculate total reads, discarded reads, and survival rate, appending results to a BioProject-specific CSV. 
#Supports single-end (se) and paired-end (pe) formats.

# ./log_trim_reads.sh $path_to_/$SRR.log format bioproject srr_accession
#------------------------------------------------------------------------------
BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"


input_file="$1"
format="$2"
bioP="$3"
sample="$4"

# Function for processing single-end reads
process_single_end() {
    total_reads_prov=$(wc -l < "$input_file")
    total_reads=$(echo "$total_reads_prov" | sed 's/ .*//')
    echo "Total reads: $total_reads"

    discarded_reads=$(grep -E ".*[[:space:]]0[[:space:]]0[[:space:]]0.*" "$input_file" | wc -l)
    echo "Discarded reads: $discarded_reads"

    survival_rate_prov=$(($total_reads - $discarded_reads))
    survival_rate_prov2=$(($survival_rate_prov * 100))
    survival_rate=$(echo "scale=5 ; $survival_rate_prov2 / $total_reads" | bc)
    echo "Survival rate: $survival_rate%"

    # Log results to CSV
    printf "%s\t%d\t%d\t%.5f\n" "$sample" "$total_reads" "$discarded_reads" "$survival_rate" >> "/$BASE_DIR/$bioP/logs/trimmomatic_log.csv"
}

# Function for processing paired-end reads
process_paired_end() {
    total_reads_prov=$(wc -l < "$input_file")
    total_reads=$(echo "$total_reads_prov" | sed 's/ .*//')
    echo "Total reads: $total_reads"

    discarded_single_reads=$(grep -E ".*[[:space:]]0[[:space:]]0[[:space:]]0.*" "$input_file" | wc -l)

    pairs=$(pcregrep -Mc '.*1[[:space:]]0[[:space:]]0[[:space:]]0[[:space:]]0\n.*2[[:space:]]0[[:space:]]0[[:space:]]0' "$input_file")

    pairs2=$((pairs * 2))
    other_pairs=$(($discarded_single_reads - $pairs2))
    other_pairs2=$((other_pairs * 2))
    discarded_reads=$((pairs2 + other_pairs2))
    echo "Discarded reads: $discarded_reads"

    survival_rate_prov=$(($total_reads - $discarded_reads))
    survival_rate_prov2=$(($survival_rate_prov * 100))
    survival_rate=$(echo "scale=5 ; $survival_rate_prov2 / $total_reads" | bc)
    echo "Survival rate: $survival_rate%"

    # Log results to CSV
    printf "%s\t%d\t%d\t%.5f\n" "$sample" "$total_reads" "$discarded_reads" "$survival_rate" >> "$BASE_DIR/$bioP/logs/trimmomatic_log.csv"
}


if [ "$format" == "se" ]; then
    process_single_end
elif [ "$format" == "pe" ]; then
    process_paired_end
fi