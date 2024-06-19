# NNMF-for-scRNAseq
Non negative matrix factorization for Oligo lineage scRNAseq data

### Analysis Pipeline for scRNA-seq Data Using NMF

This repository contains an R script (`script.R`) for analyzing single-cell RNA sequencing (scRNA-seq) data using Non-negative Matrix Factorization (NMF) with Seurat and RcppML libraries. Below is a step-by-step guide on how to reproduce the analysis:

#### Steps to Reproduce

1. **Load Required Libraries**
   - Seurat: `library(Seurat)`
   - dplyr: `library(dplyr)`
   - RcppML: `library(RcppML)`

2. **Set Seed for Reproducibility**
   - `set.seed(200)`

3. **Load Seurat Object**
   - Replace `/data/nasser/Manuscript/processedobject/ODC35_woClus8_subclust3_res0.15_NK` with your own path.
   - `ol <- readRDS("/data/nasser/Manuscript/processedobject/ODC35_woClus8_subclust3_res0.15_NK")`

4. **Set Cell Type Identities**
   - `Idents(ol) <- "CellType"`

5. **Subset the Data**
   - Choose specific cell types for analysis (`iODC`, `iOPC`, `iPPC_0`, `iPPC_1`, `iPPC_2`).
   - Uncomment and modify if subsetting by `iCEP` is required.
   - `pd <- subset(ol, idents = c("iODC", "iOPC", "iPPC_0", "iPPC_1", "iPPC_2"))`

6. **Optional: Visualize Using DimPlot**
   - Uncomment `DimPlot(pd)` to visualize the subset data.

7. **Clean the Cluster (Optional)**
   - Clean UMAP coordinates to focus on specific areas (`umap1 > -2 & umap2 > 1`).

8. **Extract Expression Matrix**
   - Extract the RNA expression matrix from the subsetted Seurat object.
   - `expression_matrix <- LayerData(object = pd, assay = "RNA", layer = "data")`
   - Remove rows with NA or null values and those with all zero values.

9. **Set Number of Clusters**
   - Determine the number of clusters based on unique cell types.
   - `num_clusters <- length(unique(pd$CellType))`

10. **Perform NMF**
    - Apply NMF to the expression matrix.
    - `nmf_result <- nmf(expression_matrix, k = num_clusters, tol = 1e-4, maxit = 500)`

11. **Extract Basis (W) and Coefficient (H) Matrices**
    - Retrieve the basis (W) and coefficient (H) matrices from the NMF result.
    - `W <- nmf_result$w`
    - `H <- nmf_result$h`

12. **Identify Most Influential Genes**
    - Determine top influential genes for each cluster.
    - Store results in `influential_genes` list.

13. **Convert Results to Data Frame**
    - Convert the list of influential genes to a tidy data frame (`influential_genes_df`).

14. **Example Visualization**
    - Generate feature plots for each cluster using influential genes.

#### Notes
- Adjust paths (`readRDS`) and parameters (`subset`, `nmf`) based on your specific data and analysis requirements.
- Ensure all necessary libraries (`Seurat`, `dplyr`, `RcppML`) are installed and loaded before running the script.

#### License
- This code is licensed under [MIT License](https://opensource.org/licenses/MIT).


