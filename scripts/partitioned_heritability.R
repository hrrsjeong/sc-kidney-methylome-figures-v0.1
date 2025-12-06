setwd("~/Projects/sciMetv2/partitioned_het/")
# Install necessary packages if you don't have them already
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

# Load the libraries
library(pheatmap)
library(RColorBrewer)


# --- 1. Configuration ---
INPUT_FILE <- "signed_heatmap_data_all.csv"
OUTPUT_FILE <- "enrichment_heatmap_annotated_tmp.pdf"
SIGNIFICANCE_THRESHOLD <- 0.05


# --- 2. Load and Prepare Data ---
# Check if the input file exists
if (!file.exists(INPUT_FILE)) {
  stop(paste("Input file not found:", INPUT_FILE))
}

# Load the matrix with signed -log10(p-values)
signed_values_matrix <- as.matrix(read.csv(INPUT_FILE, row.names = 1, check.names = FALSE))

# Clean the GWAS trait names (row names)
original_rownames <- rownames(signed_values_matrix)
cleaned_rownames <- gsub("PASS_", "", original_rownames)
cleaned_rownames <- gsub("UKB_460k.biochemistry_", "", cleaned_rownames)
rownames(signed_values_matrix) <- cleaned_rownames


# --- 3. Reconstruct P-values for Significance Testing ---
# Get the absolute value, which corresponds to the unsigned -log10(p-value)
abs_log_p_matrix <- abs(signed_values_matrix)

# Reverse the -log10 transformation to get the original p-values back
reconstructed_p_matrix <- 10^(-abs_log_p_matrix)
reconstructed_p_matrix

# --- 4. Prepare Annotations for Significant Cells ---
# Use the reconstructed p-values to determine significance
significance_matrix <- ifelse(reconstructed_p_matrix < SIGNIFICANCE_THRESHOLD, "*", "")
#border_matrix <- ifelse(reconstructed_p_matrix < SIGNIFICANCE_THRESHOLD, "black", NA)


# --- 5. Customize Heatmap Aesthetics ---
# Create a symmetric diverging color palette (Blue -> White -> Red)
palette <- rev(brewer.pal(n = 11, name = "RdBu"))

# Determine the symmetric range for the color scale
max_abs_val <- max(abs(signed_values_matrix), na.rm = TRUE)
color_breaks <- seq(-max_abs_val, max_abs_val, length.out = 101)

signed_values_matrix
# --- 6. Generate and Save the Heatmap ---d
pheatmap(
  mat               = signed_values_matrix,
  color             = colorRampPalette(palette)(100),
  breaks            = color_breaks,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  #main              = "Signed -log10(P-value) of Heritability Enrichment",
  fontsize_row      = 9,
  fontsize_col      = 9,
  display_numbers   = significance_matrix,
  fontsize_number   = 8,
  #border_color      = border_matrix,
  #filename          = OUTPUT_FILE,
  width             = 10,  # Adjust width as needed
  height            = 14   # Adjust height as needed
)
#draw()
# Print a confirmation message to the console
print(paste("✔ Publication-quality heatmap saved to:", OUTPUT_FILE))

