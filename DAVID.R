# =============================================================================
# DAVID ANALYSIS - GSE152418 (COVID-19 vs Healthy)
# =============================================================================

# =============================================================================
# STEP 1: LOAD METADATA FROM SERIES MATRIX
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 1: Loading metadata from GEO series matrix\n")
cat(rep("=", 70), "\n\n")

metadata_path <- "C:/Users/Hp/Downloads/DATA FOR EXPERIMENT/assignment_29nov/david_4d/GSE152418_series_matrix.txt"

# Read series matrix file (skip comment lines starting with !)
metadata_lines <- readLines(metadata_path)

# Find sample names
sample_line <- grep("^!Sample_title", metadata_lines, value = TRUE)
sample_names <- strsplit(sample_line, "\t")[[1]][-1]
sample_names <- gsub('"', '', sample_names)

# Find disease states
disease_line <- grep("^!Sample_characteristics_ch1.*disease state:", metadata_lines, value = TRUE)[1]
disease_states <- strsplit(disease_line, "\t")[[1]][-1]
disease_states <- gsub('"', '', disease_states)
disease_states <- gsub("disease state: ", "", disease_states)

# Create metadata dataframe
metadata <- data.frame(
  sample = sample_names,
  disease_state = disease_states,
  stringsAsFactors = FALSE
)

cat("Total samples:", nrow(metadata), "\n\n")
cat("Sample groups:\n")
print(table(metadata$disease_state))
cat("\nFirst few samples:\n")
print(head(metadata, 10))

# =============================================================================
# STEP 2: LOAD RAW COUNTS
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 2: Loading raw counts\n")
cat(rep("=", 70), "\n\n")

counts_path <- "C:/Users/Hp/Downloads/DATA FOR EXPERIMENT/assignment_29nov/david_4d/GSE152418_p20047_Study1_RawCounts.txt"
counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)

cat("Counts shape:", nrow(counts), "genes ×", ncol(counts), "samples\n")
cat("Count columns:", paste(head(colnames(counts), 5), collapse = ", "), "...\n")

# =============================================================================
# STEP 3: MATCH METADATA TO COUNTS COLUMNS
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 3: Matching metadata samples to count columns\n")
cat(rep("=", 70), "\n\n")

# Create sample to group mapping
sample_to_group <- setNames(metadata$disease_state, metadata$sample)

# Match samples - try direct match first
covid_samples <- c()
healthy_samples <- c()

for (col in colnames(counts)) {
  # Try to find matching sample in metadata
  matched <- FALSE
  for (i in 1:nrow(metadata)) {
    if (grepl(metadata$sample[i], col, fixed = TRUE) || 
        grepl(col, metadata$sample[i], fixed = TRUE)) {
      if (grepl("COVID", metadata$disease_state[i])) {
        covid_samples <- c(covid_samples, col)
      } else if (grepl("Healthy", metadata$disease_state[i])) {
        healthy_samples <- c(healthy_samples, col)
      }
      matched <- TRUE
      break
    }
  }
}

cat("COVID-19 samples:", length(covid_samples), "\n")
cat("Healthy samples:", length(healthy_samples), "\n")

if (length(covid_samples) == 0 || length(healthy_samples) == 0) {
  cat("\n❌ ERROR: Could not match samples!\n")
  cat("Count columns:", paste(head(colnames(counts)), collapse = ", "), "\n")
  cat("Metadata samples:", paste(head(metadata$sample), collapse = ", "), "\n")
  stop("Sample matching failed")
}

# =============================================================================
# STEP 4: DIFFERENTIAL EXPRESSION (COVID-19 vs Healthy)
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 4: Differential Expression - COVID-19 vs Healthy\n")
cat(rep("=", 70), "\n\n")

# Filter low counts
mean_counts <- rowMeans(counts)
counts_filt <- counts[mean_counts > 10, ]
cat("Genes after filtering (mean > 10):", nrow(counts_filt), "\n")

# Calculate fold changes
healthy_mean <- rowMeans(counts_filt[, healthy_samples]) + 1
covid_mean <- rowMeans(counts_filt[, covid_samples]) + 1
log2fc <- log2(covid_mean / healthy_mean)

# T-tests
pvals <- sapply(1:nrow(counts_filt), function(i) {
  t.test(as.numeric(counts_filt[i, covid_samples]), 
         as.numeric(counts_filt[i, healthy_samples]))$p.value
})

# Results table
results <- data.frame(
  gene_id = rownames(counts_filt),
  log2fc = log2fc,
  pvalue = pvals,
  stringsAsFactors = FALSE
)

# FDR correction
results$padj <- p.adjust(results$pvalue, method = "BH")

# Filter significant genes
sig <- results[abs(results$log2fc) > 1 & results$padj < 0.05, ]
sig <- sig[order(-abs(sig$log2fc)), ]

cat("\nSignificant genes (|log2FC| > 1, padj < 0.05):", nrow(sig), "\n")
cat("  Upregulated in COVID-19:", sum(sig$log2fc > 0), "\n")
cat("  Downregulated in COVID-19:", sum(sig$log2fc < 0), "\n")

cat("\nTop 10 upregulated genes in COVID-19:\n")
print(head(sig[sig$log2fc > 0, ], 10))

cat("\nTop 10 downregulated genes in COVID-19:\n")
print(head(sig[sig$log2fc < 0, ], 10))

# =============================================================================
# STEP 5: SAVE GENE LISTS FOR DAVID
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 5: Saving gene lists for DAVID\n")
cat(rep("=", 70), "\n\n")

# All significant genes
sig_sorted <- sig[order(-sig$log2fc), ]
write.table(sig_sorted$gene_id, 
            file = "significant_genes_for_david.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("✓ significant_genes_for_david.txt (", nrow(sig), " genes)\n", sep = "")

# Upregulated
up <- sig[sig$log2fc > 0, ]
up <- up[order(-up$log2fc), ]
write.table(up$gene_id,
            file = "upregulated_genes_for_david.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("✓ upregulated_genes_for_david.txt (", nrow(up), " genes)\n", sep = "")

# Downregulated
down <- sig[sig$log2fc < 0, ]
down <- down[order(down$log2fc), ]
write.table(down$gene_id,
            file = "downregulated_genes_for_david.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("✓ downregulated_genes_for_david.txt (", nrow(down), " genes)\n", sep = "")

# Background
write.table(rownames(counts_filt),
            file = "background_genes_for_david.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("✓ background_genes_for_david.txt (", nrow(counts_filt), " genes)\n", sep = "")

# Save full results
write.csv(results, file = "differential_expression_results.csv", row.names = FALSE)
cat("✓ differential_expression_results.csv\n")

# =============================================================================
# STEP 6: VOLCANO PLOT
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 6: Creating volcano plot\n")
cat(rep("=", 70), "\n\n")

png("volcano_plot.png", width = 10, height = 6, units = "in", res = 300)

# Base plot
plot(results$log2fc, -log10(results$pvalue),
     pch = 20, cex = 0.5, col = "gray",
     xlab = "Log2 Fold Change (COVID-19 vs Healthy)",
     ylab = "-Log10(p-value)",
     main = "Volcano Plot: COVID-19 vs Healthy")

# Significant points
sig_up <- results$log2fc > 1 & results$padj < 0.05
sig_down <- results$log2fc < -1 & results$padj < 0.05

points(results$log2fc[sig_up], -log10(results$pvalue[sig_up]),
       pch = 20, cex = 0.8, col = "red")
points(results$log2fc[sig_down], -log10(results$pvalue[sig_down]),
       pch = 20, cex = 0.8, col = "blue")

# Threshold lines
abline(h = -log10(0.05), lty = 2, col = "black")
abline(v = c(-1, 1), lty = 2, col = "black")

# Legend
legend("topright", 
       legend = c(paste("Up in COVID-19 (n=", sum(sig_up), ")", sep = ""),
                  paste("Down in COVID-19 (n=", sum(sig_down), ")", sep = "")),
       col = c("red", "blue"), pch = 20, cex = 0.8)

dev.off()
cat("✓ volcano_plot.png saved\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("ANALYSIS SUMMARY\n")
cat(rep("=", 70), "\n\n")

cat("Dataset: GSE152418 (COVID-19 vs Healthy)\n")
cat("Comparison: COVID-19 (n=", length(covid_samples), ") vs Healthy (n=", 
    length(healthy_samples), ")\n\n", sep = "")

cat("Differential Expression:\n")
cat("- Total genes analyzed:", nrow(counts_filt), "\n")
cat("- Significant genes:", nrow(sig), "(|log2FC| > 1, padj < 0.05)\n")
cat("  * Upregulated in COVID-19:", sum(sig$log2fc > 0), "\n")
cat("  * Downregulated in COVID-19:", sum(sig$log2fc < 0), "\n\n")

cat("Files generated for DAVID:\n")
cat("1. significant_genes_for_david.txt\n")
cat("2. upregulated_genes_for_david.txt\n")
cat("3. downregulated_genes_for_david.txt\n")
cat("4. background_genes_for_david.txt\n")

# =============================================================================
# NEXT STEPS
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("NEXT STEPS - UPLOAD TO DAVID\n")
cat(rep("=", 70), "\n\n")

cat("1. Go to: https://david.ncifcrf.gov/\n")
cat("2. Click 'Start Analysis'\n")
cat("3. Upload: significant_genes_for_david.txt\n")
cat("4. Select identifier: ENSEMBL_GENE_ID or OFFICIAL_GENE_SYMBOL\n")
cat("5. Select species: Homo sapiens\n")
cat("6. Click 'Submit'\n")
cat("7. Click 'Functional Annotation Tool'\n")
cat("8. Expand: GOTERM_BP_DIRECT and KEGG_PATHWAY\n")
cat("9. Take screenshots of top results\n")
cat("10. Download results\n\n")

cat("Expected enrichments for COVID-19:\n")
cat("- Immune response (GO:0006955)\n")
cat("- Inflammatory response (GO:0006954)\n")
cat("- Cytokine-mediated signaling (GO:0019221)\n")
cat("- Response to virus (GO:0009615)\n")
cat("- Type I interferon signaling\n")
cat("- Cytokine-cytokine receptor interaction (hsa04060)\n")
cat("- JAK-STAT signaling (hsa04630)\n\n")

cat(rep("=", 70), "\n\n")