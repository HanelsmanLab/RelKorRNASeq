# Load required libraries
library(edgeR)
library(stringr)
library(dplyr)

# Define the path to your files
file_path <- "Z:/projects1/Handelsman-Communities/Shruthi_Magesh/ShruthiRNASeq/htseq/"

# List all HTSeq files
file_list <- list.files(path = file_path, pattern = "*.txt", full.names = TRUE)

# Extract treatments from file names
treatments <- str_extract(basename(file_list), "f(?:-\\d+|2223)(?:-\\d+)?-(c|kora)-(t1-5|t30)")

# Debug: Print any files with missing treatments
missing_treatments <- file_list[is.na(treatments)]
if(length(missing_treatments) > 0) {
  print("Files with missing treatments:")
  print(missing_treatments)
}

# Clean treatments: Remove the first number following "f-" or "f2223-"
treatments_clean <- gsub("f-\\d+-", "f-", treatments)
treatments_clean <- gsub("f2223-\\d+-", "f2223-", treatments_clean)

# Create a factor for treatments
group <- factor(treatments_clean, levels = unique(treatments_clean))

# Read each file and store data in a list
data_list <- lapply(seq_along(file_list), function(i) {
  data <- read.table(file_list[i], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("Gene", paste0("Sample", i))
  return(data)
})

# Combine all data frames by "Gene" column
combined_data <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), data_list)

# Rename columns with cleaned treatment names
colnames(combined_data) <- c("Gene", treatments_clean)

# Convert to matrix
count_matrix <- as.matrix(combined_data[, -1])
rownames(count_matrix) <- combined_data$Gene

# Remove rows with special categories
special_categories <- c("__no_feature", "__not_aligned", "__too_low_aQual", "__ambiguous", "__alignment_not_unique")
count_matrix <- count_matrix[!rownames(count_matrix) %in% special_categories, ]

# Generate dynamic column names for each treatment
num_samples_per_treatment <- table(treatments_clean)
column_names <- unlist(lapply(unique(treatments_clean), function(treatment) {
  paste0(treatment, ".", seq_len(num_samples_per_treatment[treatment]))
}))

print(column_names)

# Ensure column names match the count matrix
if (length(column_names) != ncol(count_matrix)) {
  stop("Mismatch between generated column names and number of columns in the count matrix")
}

# Assign column names to the count matrix
colnames(count_matrix) <- column_names

# Create DGEList object and design matrix
dge <- DGEList(counts = count_matrix, group = group)
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("-", "_", levels(group))  # Replace hyphens with underscores

# Normalize data and filter lowly expressed genes
dge <- calcNormFactors(dge)
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the model
fit <- glmQLFit(dge, design)

# Define contrasts for the analysis
contrasts <- makeContrasts(
  "Line51f_c_t30_vs_f_c_t1_5" = "f_c_t30 - f_c_t1_5", #A
  "Line51f_kora_t30_vs_f_kora_t1_5" = "f_kora_t30 - f_kora_t1_5", #D
  "Line51f_kora_t30_vs_f_c_t30" = "f_kora_t30 - f_c_t30", #C
  "Line51f_kora_t1_5_vs_f_c_t1_5" = "f_kora_t1_5 - f_c_t1_5", #B
  "Line51f2223_c_t30_vs_f2223_c_t1_5" = "f2223_c_t30 - f2223_c_t1_5", #I
  "Line51f2223_c_t30_vs_f_c_t30" = "f2223_c_t30 - f_c_t30", #F
  "Line51f2223_kora_t30_vs_f_kora_t30" = "f2223_kora_t30 - f_kora_t30",  #H
  "Line51f2223_kora_t1_5_vs_f2223_c_t1_5" = "f2223_kora_t1_5 - f2223_c_t1_5",  #J
  "Line51f2223_kora_t1_5_vs_f_kora_t1_5" = "f2223_kora_t1_5 - f_kora_t1_5", #G
  "Line51f2223_kora_t30_vs_f2223_c_t30" = "f2223_kora_t30 - f2223_c_t30", #K
  "Line51f2223_kora_t30_vs_f2223_kora_t1_5" = "f2223_kora_t30 - f2223_kora_t1_5", #L
  "Line51f2223_c_t1_5_vs_f_c_t1_5" = "f2223_c_t1_5 - f_c_t1_5", #E
  levels = design
)


# Initialize results_df_list as an empty list
results_df_list <- list()

# Apply each contrast and store results with proper names
for (contrast_name in colnames(contrasts)) {
  # Perform GLM QL test for the current contrast
  result <- glmQLFTest(fit, contrast = contrasts[, contrast_name])
  
  # Extract result table
  result_df <- topTags(result, n = Inf)$table
  
  # Assign the result to the list with the contrast name as the key
  results_df_list[[contrast_name]] <- result_df
}

# Check if the results_df_list has been populated correctly
print(names(results_df_list))  # This should now print the contrast names

# Save results with meaningful filenames based on the contrast name
for (name in names(results_df_list)) {
  file_name <- paste0(name, "_results.csv")
  
  # Debugging: Print filename and check list content
  print(paste("Attempting to save results to:", file_name))
  if (nrow(results_df_list[[name]]) > 0) {
    write.csv(results_df_list[[name]], file = file_name, row.names = TRUE)
    print(paste("Successfully saved:", file_name))
  } else {
    print(paste("No data to save for:", name))
  }
}

# Verify files in directory
print(list.files())
