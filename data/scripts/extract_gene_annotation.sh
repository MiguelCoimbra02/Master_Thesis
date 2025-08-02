#!/bin/bash
# Usage: ./extract_gene_annotation.sh gene_list.txt annotation.gtf output.csv
# or
# Usage: ./extract_gene_annotation.sh -> use default values

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
GTF_DIR="$BASE_DIR/genome/genome_annotation"


mkdir -p "$GTF_DIR"

# Check if the GTF directory is empty
if [ -z "$(ls -A "$GTF_DIR")" ]; then
    echo "The GTF directory is empty. Run this command if you haven't yet: ./star_fc.sh BioProject format(se or pe)"
    exit 1
fi

GENE_LIST="${1:-$BASE_DIR/seidr_output/genes.txt}"
GTF_FILE="${2:-$GTF_DIR/cork_oak_genome_annotation.gtf}"
OUTPUT_FILE="${3:-$GTF_DIR/gene_functional_annotation.csv}"

# Create or clear the output file and add a header
echo -e "GeneID\tDescription" > "$OUTPUT_FILE"

# Process each gene and extract its description
while IFS= read -r GENE; do
    # Find the annotation line with the specific gene ID
    ANNOTATION=$(grep -P "\tgene\t.*gene_id \"$GENE\";" "$GTF_FILE")
    
    if [ -n "$ANNOTATION" ]; then
        # Extract only the description
        DESCRIPTION=$(echo "$ANNOTATION" | grep -o 'description "[^"]*"' | cut -d'"' -f2)
        
        # Check if the description is "uncharacterized" or similar and replace with a default message
        if [[ "$DESCRIPTION" =~ uncharacterized ]]; then
            DESCRIPTION="No Description available"
        fi
        
        # Provide a default message if no description is found
        [ -z "$DESCRIPTION" ] && DESCRIPTION="No description found"
    else
        DESCRIPTION="No annotation found"
    fi
    
    # Append the gene and description to the output file
    echo -e "$GENE\t$DESCRIPTION" >> "$OUTPUT_FILE"
done < "$GENE_LIST"

echo "Gene annotations have been saved to $OUTPUT_FILE"
