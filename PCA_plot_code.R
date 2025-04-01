# Load required libraries
library(ggplot2)
library(plotly)
library(dplyr)

data <- read.csv("shruthi_cpm_valuesNEWCODE.csv", row.names = 1)  # Ensure row names are set properly

# Transpose the dataset so that samples are rows and genes are columns
transposed_data <- as.data.frame(t(data))  # Convert matrix to data frame

# Convert all columns to numeric (if transposing introduced character data)
transposed_data[] <- lapply(transposed_data, function(x) as.numeric(as.character(x)))

# Remove genes (columns) with zero variance
transposed_data <- transposed_data[, apply(transposed_data, 2, var, na.rm = TRUE) > 0]

# Remove rows with missing values
transposed_data <- na.omit(transposed_data)

# Perform PCA on the filtered data
pca_result <- prcomp(transposed_data, center = TRUE, scale. = TRUE)

# Extract PCA scores for samples
pca_data <- as.data.frame(pca_result$x)

# Add sample names as a new column
pca_data$Sample <- rownames(pca_data)

# Extract general treatment names (removing biological replicate numbers)
pca_data$Treatment <- gsub("^(f2223|f)\\.\\d+", "\\1", pca_data$Sample)  # Remove biological rep number
pca_data$Treatment <- gsub("_Fj$", "", pca_data$Treatment)  # Remove "_Fj" at the end


# Check extracted treatments
print(unique(pca_data$Treatment))  # Debugging step

# Generate colors for each unique treatment
unique_treatments <- unique(pca_data$Treatment)
color_palette <- setNames(scales::hue_pal()(length(unique_treatments)), unique_treatments)

# Check if colors are assigned correctly
print(color_palette)  # Debugging step

# Create PCA plot with treatment groups and clustering ellipses
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3) +  # Scatter points
  xlab(paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)")) +
  scale_color_manual(values = color_palette) +  # Assign fixed colors
  theme_minimal() +
  ggtitle("PCA Plot of Treatments with Clustering")

# Convert to interactive plot
ggplotly(pca_plot)

write.csv(pca_data, "PCA_scores.csv", row.names = FALSE)


