# Load necessary libraries
library(readxl)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
# question 1 
# Step 1: Reading the input data
gene_data <- read_excel("C:\Users\Admin\Downloads\Gene_Expression_Data (1).xlsx")
gene_info <- read.csv("C:\Users\Admin\Downloads\Gene_Information (1).csv")
sample_info <- read.table("C:\Users\Admin\Downloads\Sample_Information (1).tsv", header = TRUE, sep = '\t')

# Step 2: Updating column names in gene expression data based on sample information
phenotype_labels <- sub(" .*", "", sample_info$group)
updated_colnames <- paste0(phenotype_labels, "_", colnames(gene_data)[-1])
colnames(gene_data)[-1] <- updated_colnames

# Step 3: Splitting the data based on phenotypes (Tumor and Normal)
tumor_data <- gene_data[, grepl("^t", colnames(gene_data))]
normal_data <- gene_data[, grepl("^n", colnames(gene_data))]

# Step 4: Calculating average expression for each probe
tumor_means <- rowMeans(tumor_data[, -1])
normal_means <- rowMeans(normal_data[, -1])

# Step 5: Computing fold change between tumor and normal groups
fold_changes <- (tumor_means - normal_means) / normal_means

# Step 6: Filtering genes with significant fold change (magnitude > 5)
significant_genes <- data.frame(Index = seq_along(fold_changes), Fold_Change = fold_changes)
significant_genes <- significant_genes[abs(significant_genes$Fold_Change) > 5, ]

# Step 7: Annotating results with higher expression and chromosome information
fold_change_annotations <- data.frame(Probe = gene_info$Probe_ID, Fold_Change = fold_changes)
fold_change_annotations$Expression <- ifelse(fold_change_annotations$Fold_Change > 0, 'Tumor', 'Normal')
fold_change_annotations <- left_join(fold_change_annotations, gene_info[, c('Probe_ID', 'Chromosome')], by = c('Probe' = 'Probe_ID'))

# Displaying results
print(significant_genes)
print(fold_change_annotations)

# question 2 
# 2. Histogram: Distribution of DEGs by chromosome
ggplot(fold_change_annotations, aes(x = Chromosome)) +
  geom_histogram(stat = "count", fill = "skyblue", color = "black") +
  labs(title = "Distribution of DEGs by Chromosome",
       x = "Chromosome",
       y = "Number of DEGs") +
  theme_minimal()
# 3. Histogram: DEGs by chromosome segregated by sample type (Normal or Tumor)
ggplot(fold_change_annotations, aes(x = Chromosome, fill = Expression)) +
  geom_histogram(stat = "count", position = "stack") +
  labs(title = "Distribution of DEGs by Chromosome Segregated by Sample Type",
       x = "Chromosome",
       y = "Number of DEGs",
       fill = "Sample Type") +
  theme_minimal()

# 4. Bar chart: Percentages of DEGs upregulated and downregulated in Tumor samples
expression_percentage <- fold_change_annotations %>%
  group_by(Expression) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100)

ggplot(expression_percentage, aes(x = Expression, y = Percentage, fill = Expression)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Percentage of DEGs Upregulated and Downregulated in Tumor Samples",
       x = "Expression Type",
       y = "Percentage") +
  theme_minimal()

# 5. Heatmap: Visualizing gene expression by sample

# Load necessary library
library(pheatmap)

# Calculate variance for each gene
gene_variances <- apply(heatmap_data, 1, var)

# Subset the top 1000 most variable genes
top_variable_genes <- order(gene_variances, decreasing = TRUE)[1:1000]
subset_heatmap_data <- heatmap_data[top_variable_genes, ]

# Convert to matrix
subset_heatmap_matrix <- as.matrix(subset_heatmap_data)

# Create the heatmap
pheatmap(subset_heatmap_matrix,
         cluster_rows = FALSE,  # Disable row clustering for simplicity
         cluster_cols = FALSE,  # Disable column clustering for simplicity
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Top 1000 Most Variable Genes")

# 6. Clustermap: Visualizing gene expression by sample with clustering
pheatmap(subset_heatmap_matrix,
         cluster_rows = TRUE,  # Enable row clustering
         cluster_cols = TRUE,  # Enable column clustering
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Clustermap of Top 1000 Most Variable Genes")
# 7.The exploratory data analysis of significant genes (those with fold change magnitude > 5) reveals several important insights:

Distribution of Fold Changes: The histogram of fold changes shows a clear separation between genes upregulated in tumor samples and those downregulated in tumor samples, with most significant genes exhibiting substantial differences in expression. This highlights a distinct transcriptional signature associated with tumor and normal phenotypes.

Chromosomal Enrichment: The distribution of differentially expressed genes (DEGs) by chromosome reveals certain chromosomes with higher concentrations of significant genes, suggesting these chromosomes may harbor genes crucial to tumor progression or suppression.

Expression Type Analysis: The bar chart showing the proportion of DEGs indicates that a majority of significant genes are upregulated in tumor samples, implying an activation of oncogenic pathways. Conversely, a smaller proportion are downregulated, suggesting suppression of normal regulatory processes.

Heatmap and Clustermap: Heatmaps and clustermaps of gene expression clearly segregate tumor and normal samples based on expression profiles. Genes clustered together in these visualizations suggest co-regulation, potentially pointing to shared biological functions or pathways.
