# scripts/rna_analysis.R

# Load necessary libraries
# Bioconductor packages
library(tximport)
library(DESeq2)
library(PCAtools)
library(GenomicRanges)

# CRAN packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(plotly)
library(ape) # For phylogenetic trees, if you want that
library(yaml) # To read the config.yaml file

# ==============================================================================
# 1. READ CONFIGURATION AND SETUP DATA
# ==============================================================================
# Read the Snakefile config to get sample information
config <- read_yaml("config.yaml")
samples <- names(config$samples)
files <- file.path("results/counts", paste0(samples, "_kallisto.tsv"))
names(files) <- samples

# Create a sample table
coldata <- data.frame(
    row.names = samples,
    condition = c("ctrl", "ctrl", "exp", "exp") # Assumes the order in config.yaml
)

# Use tximport to read kallisto output
txi <- tximport(files, type = "kallisto", txOut = TRUE)

# ==============================================================================
# 2. DIFFERENTIAL EXPRESSION ANALYSIS (DESeq2)
# ==============================================================================
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)
dds <- DESeq(dds)
res <- results(dds)

# Save the results
write.csv(as.data.frame(res), "results/R_analysis/deseq2_results.csv")

# ==============================================================================
# 3. PLOTTING
# ==============================================================================
# 3.1. PCA Plot
# First, get the variance-stabilized data
vst <- vst(dds, blind = FALSE)
pca_data <- pca(assay(vst), metadata = coldata, removeVar = 0.1)

# Plot using PCAtools
png("results/R_analysis/pca_plot.png")
p <- biplot(pca_data, x = 'PC1', y = 'PC2', lab = coldata$condition)
print(p)
dev.off()

# 3.2. MA Plot (Mean-Average)
png("results/R_analysis/ma_plot.png")
plotMA(res, main = "MA plot", ylim = c(-2, 2))
dev.off()

# 3.3. Volcano Plot
res_df <- as.data.frame(res) %>%
    na.omit() %>%
    mutate(
        gene_name = rownames(.)
    ) %>%
    mutate(
        significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant")
    )

p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant, text = gene_name)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    labs(title = "Volcano Plot", x = "log2(Fold Change)", y = "-log10(Adjusted p-value)") +
    theme_minimal()

png("results/R_analysis/volcano_plot.png")
print(p)
dev.off()
# Use plotly for an interactive version:
# plotly_plot <- ggplotly(p)
# htmlwidgets::saveWidget(plotly_plot, "results/R_analysis/volcano_plot.html")


# 3.4. P-value Distribution Plot
png("results/R_analysis/p_value_plot.png")
hist(res$pvalue[res$baseMean > 1], breaks = 20, col = "skyblue", border = "black",
     main = "P-value Distribution", xlab = "p-value")
dev.off()


# 3.5. Gene Regulation Plot (Top 20 differentially expressed genes)
# Assuming you want to plot the expression of the top 20 DEGs.
top_genes <- res_df %>%
    arrange(padj) %>%
    head(20) %>%
    rownames()

normalized_counts <- counts(dds, normalized = TRUE)
top_counts <- normalized_counts[top_genes, ]
top_counts_df <- as.data.frame(top_counts) %>%
    rownames_to_column(var = "gene") %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
    left_join(coldata %>% rownames_to_column(var = "sample"), by = "sample")

p_gene_reg <- ggplot(top_counts_df, aes(x = gene, y = count, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Top 20 Differentially Expressed Contigs", x = "Contig ID", y = "Normalized Count") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

png("results/R_analysis/gene_regulation_plot.png")
print(p_gene_reg)
dev.off()

# 3.6. Contig Heatmap of Correlation
# Heatmap of the top 50 DEGs
top_50_genes <- res_df %>%
    arrange(padj) %>%
    head(50) %>%
    rownames()

heatmap_data <- assay(vst)[top_50_genes, ]
pheatmap_plot <- pheatmap(heatmap_data, scale = "row",
                          show_rownames = FALSE,
                          annotation_col = coldata,
                          main = "Heatmap of Top 50 DE Contigs")

png("results/R_analysis/heatmap.png")
print(pheatmap_plot)
dev.off()
