#!/usr/bin/Rscript
#This is a standalone R file since i was having problems with some commands when done with relative paths
###
# Get the path to the script directory

#Change ths directory for the wanted directory.
script_dir <- "/home/mlcoimbraj"

# Define relative paths based on the script directory
data_dir <- file.path(script_dir, "data")
bioP_file <- file.path(data_dir, "bioProjects_info", "COak_BioP_list.txt")

# Check if the required directories exist
if (!dir.exists(data_dir)) {
  stop("Data directory not found: ", data_dir)
}

if (!dir.exists(file.path(data_dir, "bioProjects_info"))) {
  stop("bioProjects_info directory not found.")
}

# Read the content of COak_BioP_list.txt
if (file.exists(bioP_file)) {
  bioP_list <- readLines(bioP_file) # Read the list of bioProjects
} else {
  stop("BioProjects file not found.")
}

# List all files in the directory that start with "PRJ"
prj_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
prj_files <- prj_dirs[basename(prj_dirs) %in% bioP_list]

# Initialize the Raw_Counts table
Raw_Counts <- data.frame()

# Loop through each file
for (prj_dir in prj_files) {
  prj_name <- basename(prj_dir)
  raw_counts_file <- file.path(prj_dir, "logs", paste0(prj_name, "_raw_counts.txt"))
  
  if (file.exists(raw_counts_file)) {
    # Read and clean the data
    table_prj_name <- read.delim(file = raw_counts_file, comment.char = "#", row.names = 1)
    table1_prj_name_clean <- table_prj_name[-c(1:5)]
    
    # Clean column names
    names(table1_prj_name_clean) <- gsub(x = names(table1_prj_name_clean), pattern = "_sorted.bam", replacement = "")
    names(table1_prj_name_clean) <- gsub(x = names(table1_prj_name_clean), pattern = ".*\\.", replacement = "")
    
    
    # Merge with the Raw_Counts table
    if (nrow(Raw_Counts) == 0) {
      Raw_Counts <- table1_prj_name_clean
    } else {
      Raw_Counts <- merge(Raw_Counts, table1_prj_name_clean, by = "row.names", all = TRUE)
      row.names(Raw_Counts) <- Raw_Counts$Row.names
      Raw_Counts$Row.names <- NULL
    }
  }
}
write.table(Raw_Counts, file = file.path(data_dir, "Raw_Counts.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

cat("Raw_Counts table has been saved to Raw_Counts.txt\n")
