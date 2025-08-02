#!/bin/bash
#Script used to extract promoter regions (2kb - default value) of all genes within the GCN
#Necessary inputs: $1 - size_of_promoter $2 - annotation_file.gff $3 - genome_file.fna 
#Run command: ./get_all_promoters_sequences.sh $1 $2 $3 or ./get_all_promoters_sequences.sh

# GTF file format (Gene Transfer Format)
# Columns:
# 1. Sequence (Chromosome/Scaffold) - e.g., NW_027069468.1
# 2. Source (Annotation Source) - e.g., Gnomon (NCBI gene prediction tool)
# 3. Feature Type - e.g., gene, transcript, exon, CDS
# 4. Start Position - Start coordinate of the feature (1-based)
# 5. End Position - End coordinate of the feature (inclusive)
# 6. Score - Usually ".", meaning no score available
# 7. Strand - "+" (forward) or "-" (reverse) strand
# 8. Frame - 0, 1, 2 (for CDS) or "." for other features
# 9. Attributes - Key-value pairs with gene/transcript IDs, annotations, etc.

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p $BASE_DIR/cis_elements/genes
mkdir -p $BASE_DIR/cis_elements/scaffolds
mkdir -p $BASE_DIR/cis_elements/promoters
mkdir -p $BASE_DIR/cis_elements/TF_motifs/fasta

GENES_DIR=$BASE_DIR/cis_elements/genes
SCAFFOLD_DIR=$BASE_DIR/cis_elements/scaffolds
PROMOTOR_DIR=$BASE_DIR/cis_elements/promoters
# Set default values
DEFAULT_PROMOTER_SIZE=2000
DEFAULT_ANNOTATION_GFF=$BASE_DIR/genome/genome_annotation/cork_oak_genome_annotation.gff
DEFAULT_GENOME=$BASE_DIR/genome/cork_oak_genome.fna


# Assign user input or default values
PROMOTER_SIZE=${1:-$DEFAULT_PROMOTER_SIZE}
ANNOTATION_FILE_GFF=${2:-$DEFAULT_ANNOTATION_GFF}
GENOME_FILE=${3:-$DEFAULT_GENOME}

echo "Using promoter size: $PROMOTER_SIZE pb"
echo "Using annotation file: $ANNOTATION_FILE_GFF"
echo "Using genome file: $GENOME_FILE"


#Getting a bed file with all cork oak gene coordinates (start to finish, excluding trnas and mrnas)


grep -P '\tgene\t' $ANNOTATION_FILE_GFF >  $GENES_DIR/genes_seqs.gff
awk '{print $1,$4,$5,$9,$7}' $GENES_DIR/genes_seqs.gff > $GENES_DIR/genes_seqs_1.gff
sed 's/ID=gene-//g' $GENES_DIR/genes_seqs_1.gff > $GENES_DIR/genes_seqs_2.gff
sed 's/;.* / /g' $GENES_DIR/genes_seqs_2.gff > $GENES_DIR/genes_seqs_3.bed
sed 's/ /\t/g' $GENES_DIR/genes_seqs_3.bed > $GENES_DIR/genes_seqs.bed

#rm $GENES_DIR/genes_seqs.gff
#rm $GENES_DIR/genes_seqs_2.gff
#rm $GENES_DIR/genes_seqs_3.bed

#Getting chromosome/scafold sizes in a fasta file
grep -P '\tregion\t' $ANNOTATION_FILE_GFF > $SCAFFOLD_DIR/scaffolds_size.gff
awk '{print $1,$5}' $SCAFFOLD_DIR/scaffolds_size.gff > $SCAFFOLD_DIR/scaffolds.fa
sed 's/ /\t/g' $SCAFFOLD_DIR/scaffolds.fa > $SCAFFOLD_DIR/scaffolds_size.fa

#rm scaffolds_size.gff
#rm $SCAFFOLD_DIR/scaffolds.fa
#rm $SCAFFOLD_DIR/scaffold_size_3.gtf

#Getting promoter region coordinates of all genes with "$1" kb (default=2kb)
grep '+'  $GENES_DIR/genes_seqs.bed > $PROMOTOR_DIR/CorkOak_genes_p.bed
bedtools flank -i $PROMOTOR_DIR/CorkOak_genes_p.bed -g $SCAFFOLD_DIR/scaffolds_size.fa -l $PROMOTER_SIZE -r 0 > $PROMOTOR_DIR/promotor_"$PROMOTER_SIZE"_pb.bed
echo "Flank of + reading written to promotor_"$PROMOTER_SIZE"_pb.bed"

grep '-'  $GENES_DIR/genes_seqs.bed > $PROMOTOR_DIR/CorkOak_genes_n.bed
bedtools flank -i $PROMOTOR_DIR/CorkOak_genes_n.bed -g $SCAFFOLD_DIR/scaffolds_size.fa -l 0 -r $PROMOTER_SIZE >> $PROMOTOR_DIR/promotor_"$PROMOTER_SIZE"_pb.bed
echo "Flank of - reading written to promotor_"$PROMOTER_SIZE"_pb.bed"

sed 's/ .*//g' $GENOME_FILE > $SCAFFOLD_DIR/scaffolds.fna
bedtools getfasta -fi $SCAFFOLD_DIR/scaffolds.fna -bed $PROMOTOR_DIR/promotor_"$PROMOTER_SIZE"_pb.bed -fo $PROMOTOR_DIR/promotor_"$PROMOTER_SIZE"_sequences.bed

echo "All files were created!"
