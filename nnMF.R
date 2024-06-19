# Load necessary libraries
library(Seurat)
library(dplyr)
library(RcppML)

# Set seed for reproducibility
set.seed(200)

# Load Seurat object
ol <- readRDS("/data/nasser/Manuscript/processedobject/ODC35_woClus8_subclust3_res0.15_NK")

# Set cell type identities
Idents(ol) <- "CellType"

# Subset the data
pd <- subset(ol, idents = c("iODC", "iOPC", "iPPC_0", "iPPC_1", "iPPC_2"))


# Clean the cluster (optional)
umap <- as.data.frame(Embeddings(object = pd[["umap"]]))
colnames(umap) <- c("UMAP_1", "UMAP_2")
pd$umap1 <- umap$UMAP_1
pd$umap2 <- umap$UMAP_2
pd <- subset(x = pd, subset = umap1 > -2 & umap2 > 1)

# Extract the expression matrix
expression_matrix <- LayerData(object = pd, assay = "RNA", layer = "data")
expression_matrix <- as.matrix(expression_matrix)

# Remove rows with NA or null values
expression_matrix <- expression_matrix[rowSums(is.na(expression_matrix)) == 0, ]
expression_matrix <- expression_matrix[rowSums(expression_matrix == 0) != ncol(expression_matrix), ]

# Set the number of clusters/components
num_clusters <- length(unique(pd$CellType))

# Perform NMF
nmf_result <- nmf(expression_matrix, k = num_clusters, tol = 1e-4, maxit = 500)

# Extract basis (W) and coefficient (H) matrices
W <- nmf_result$w
H <- nmf_result$h

# Identify most influential genes
influential_genes <- list()

for (i in 1:num_clusters) {
  gene_weights <- W[, i]
  top_genes <- order(gene_weights, decreasing = TRUE)[1:10]  # Top 10 genes
  influential_genes[[paste("CellType", i)]] <- rownames(expression_matrix)[top_genes]
}

# Convert the list to a data frame
influential_genes_df <- do.call(rbind, lapply(names(influential_genes), function(x) {
  data.frame(Cluster = x, Gene = influential_genes[[x]])
}))

# Example visualization for each cluster
for (i in 1:num_clusters) {
  cluster_name <- paste("CellType", i)
  FeaturePlot(pd, features = influential_genes[[cluster_name]], title = cluster_name)
}
