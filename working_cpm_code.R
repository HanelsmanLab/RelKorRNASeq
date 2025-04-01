# Load required libraries
library(edgeR)

# Read the count matrix, assuming first row is treatments and first column is gene IDs
file_paths <- list.files("Z:/projects1/Handelsman-Communities/Shruthi_Magesh/ShruthiRNASeq/htseq/", full.names = TRUE)

# Initialize a list to store count data
count_list <- list()

# Initialize a vector for gene names
gene_names <- NULL

# Loop through files and read count data
for (file in file_paths) {
  # Read HTSeq count data (assuming the first column is GeneID, second is count)
  counts <- read.table(file, header = FALSE, row.names = 1)
  
  # Store the gene names from the first column if it's the first file
  if (is.null(gene_names)) {
    gene_names <- rownames(counts)
  }
  
  # Add count data (second column) to the list
  count_list[[file]] <- counts[, 1]  # Assuming count is in the second column
}

# Combine all counts into a single matrix
count_matrix <- do.call(cbind, count_list)

# Assign gene names as row names for the count matrix
rownames(count_matrix) <- gene_names

# Get sample names from the file names (treatment information)
sample_names <- gsub(".*/(.*).txt", "\\1", file_paths)  # Adjust pattern based on your file naming convention

# Assign column names (sample names) to the count matrix
colnames(count_matrix) <- sample_names

# Create the group factor based on treatment information
group <- factor(sample_names)  # Treatments encoded in filenames

# Create a DGEList object
dge <- DGEList(counts = count_matrix, group = group)

# Calculate CPM (Counts Per Million)
cpm_values <- cpm(dge)

# Convert CPM values to a data frame
cpm_df <- as.data.frame(cpm_values)

# Add gene names as the first column
cpm_df$Gene <- rownames(cpm_values)

# Reorder columns so "Gene" comes first
cpm_df <- cpm_df[, c("Gene", colnames(cpm_df)[-ncol(cpm_df)])]

# Save the CPM values to a CSV file
write.csv(cpm_df, file = "Z:/projects1/Handelsman-Communities/Austin_Hall/2024TransATRNASeq/shruthi_cpm_valuesNEWCODE.csv", row.names = FALSE)

# Print some of the CPM values to verify
head(cpm_df)
