#!/bin/bash

# Set BASE_DIR to 'data' directory relative to the script's location
BASE_DIR="$(cd "$(dirname "$0")/../data" && pwd)"

# Base directories relative to BASE_DIR
base_dirs=(
  "$BASE_DIR/BioProjects_info"
  "$BASE_DIR/cis_elemets/fimo_out"
  "$BASE_DIR/cis_elemets/targeted_Net"
  "$BASE_DIR/cis_elemets/genes"
  "$BASE_DIR/cis_elemets/promoters"
  "$BASE_DIR/cis_elemets/scaffolds"
  "$BASE_DIR/cis_elemets/TF_motifs/fasta"
  "$BASE_DIR/connectTF/output"
  "$BASE_DIR/connectTF/connectTF_Targeted"
  "$BASE_DIR/connectTF/connectTF_Targeted/output_targeted"
  "$BASE_DIR/dap_seq/data"
  "$BASE_DIR/dap_seq/data_targeted"
  "$BASE_DIR/genome/genome_annotation"
  "$BASE_DIR/genome/genome_index"
  "$BASE_DIR/proteome"
  "$BASE_DIR/scripts"
  "$BASE_DIR/seidr_output"
)
#REVER QUANDO INSERIR DIRS DA TOOL


# Create directories relative to BASE_DIR
for dir in "${base_dirs[@]}"; do
  if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
    echo "Created: $dir"
  else
    echo "Already exists: $dir"
  fi
done

echo "Directory structure setup completed."
