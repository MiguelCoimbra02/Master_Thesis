#!/bin/bash

# Set DATA_DIR to 'data' directory relative to the script's location
DATA_DIR="$(cd "$(dirname "$0")/../data" && pwd)"
DOCKER_DIR="$(cd "$(dirname "$0")/../../viz_tool/docker_dir" && pwd)"


# Base directories relative to DATA_DIR
dirs=(
  "$DATA_DIR/BioProjects_info"
  "$DATA_DIR/cis_elements/fimo_out"
  "$DATA_DIR/cis_elements/targeted_Net"
  "$DATA_DIR/cis_elements/genes"
  "$DATA_DIR/cis_elements/promoters"
  "$DATA_DIR/cis_elements/scaffolds"
  "$DATA_DIR/cis_elements/TF_motifs/fasta"
  "$DATA_DIR/connectTF/output"
  "$DATA_DIR/connectTF/connectTF_Targeted"
  "$DATA_DIR/connectTF/connectTF_Targeted/output_targeted"
  "$DATA_DIR/dap_seq/data"
  "$DATA_DIR/dap_seq/data_targeted"
  "$DATA_DIR/genome/genome_annotation"
  "$DATA_DIR/genome/genome_index"
  "$DATA_DIR/proteome"
  "$DATA_DIR/scripts"
  "$DATA_DIR/seidr_output"
  "$DOCKER_DIR/data"
  "$DOCKER_DIR/static/images"
  "$DOCKER_DIR/templates"
)
# Create directories relative to dirs
for dir in "${dirs[@]}"; do
  if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
    echo "Created: $dir"
  else
    echo "Already exists: $dir"
  fi
done

echo "Directory structure setup completed."
