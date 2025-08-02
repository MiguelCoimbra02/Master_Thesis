#!/usr/bin/Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))

#install.packages("argparse")
#install.packages("RColorBrewer")

#install.packages("BiocManager")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("DESeq2")

library("RColorBrewer")
library("argparse")
#suppressPackageStartupMessages(library("DESeq2"))

parser <- ArgumentParser(description="Get the expression.tsv and genes.txt from Raw_Counts.txt file ")
parser$add_argument("-f", "--rawcountsfile", help="Raw_Counts.txt table delimited by spaces")
parser$add_argument("-o", "--output", metavar="STR", help="Output files (default = seidr_output)")

args <- parser$parse_args()
rawcountsfile <- args$rawcountsfile
output <- args$output

suppressPackageStartupMessages(library("DESeq2"))
Raw_Counts <- read.delim(file = rawcountsfile, comment.char = "#", row.names = 1)

#Since i dont have metadata object i need to supply DESeq2 a Data Frame with any content
dummy_meta <- data.frame(N = seq_along(Raw_Counts)) # will generate an integer per row of Raw_Counts
#DESeq2 perform the differential expression analysis without any experimental design or grouping, since my 'metadata' soes not exist
dds <- DESeqDataSetFromMatrix(Raw_Counts, dummy_meta, ~1)

#removing genes with a mean less than 10 for each sample
keep <- rowMeans(counts(dds)) >= 10
#only the genes for wich keep is TRUE will be kept
dds <- dds[keep,]

#apply variance stabilization to the data in the dds (DESeqDataSet) object
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
head(vsd)
#change genes to columns and samples to rows
vst <- t(assay(vsd))

#REMOVING GENES WITH NO VARIATION ACROSS THE SAMPLES
vars <- apply(vst, 2, var)
filt_id <- which(is.finite(vars))
vst <- vst[, filt_id]

#normalization centering the median
medians <- apply(vst, 1, median)
vst <- sweep(vst, MARGIN=1, FUN='-', STATS=medians)

#distance matrix heatmap
sampleDists <- dist(vst)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#######

# Get the path to the current working directory (or the script directory)
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

# Define relative paths for output files based on script location
expression_file <- file.path(script_dir, "data", output, "expression.tsv")
genes_file <- file.path(script_dir, "data", output, "genes.txt")

#write.matrix is similar to write.table but faster for big matrices. unname removes headers
#file with column headers; the genes themselves
MASS::write.matrix(x = unname(vst), sep = '\t', file = expression_file)
write(colnames(vst), file = genes_file)