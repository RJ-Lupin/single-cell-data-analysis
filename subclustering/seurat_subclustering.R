# Load the required packages
library(Seurat)
library(dplyr)
library(ggplot2)

# @param seurat_object: A Seurat object containing the single-cell RNA-seq data
# @param seed: An integer for setting the random seed (default: 12345)
# @param n_variable_features: Number of variable features to select (default: 2000)
# @param pca_dims: Number of principal components to compute and use (default: 20)
# @param idents_to_keep: Identities to maintain or remove for subsetting (default: NULL)
# @param invert: Logical, whether to invert the selection of identities (default: FALSE)
# @param ident_column: Metadata column name to set as the active identity class (default: NULL)

seurat_subclustering <- function(seurat_object, 
                                 seed = 12345,
                                 n_variable_features = 2000,
                                 pca_dims = 50,
                                 idents_to_keep = NULL,
                                 invert = FALSE,
                                 ident_column = NULL) {
  
  set.seed(seed)
  
  # Set the active identity class from a metadata column if specified
  if (!is.null(ident_column)) {
    Idents(seurat_object) <- ident_column
  }
  
  # Subset the Seurat object based on specified identities
  if (!is.null(idents_to_keep)) {
    seurat_object <- subset(seurat_object, idents = idents_to_keep, invert = invert)
  }
  
  # Standard Seurat workflow
  seurat.subset <- seurat_object %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = n_variable_features) %>%
    ScaleData() %>%
    RunPCA()
  
  # Elbow plot for PCA
  print(ElbowPlot(seurat.subset, ndims = 50))
  
  # Downstream analysis with PCA results
  seurat.subset <- seurat.subset %>%
    RunUMAP(dims = 1:pca_dims) %>%
    FindNeighbors(dims = 1:pca_dims) %>%
    FindClusters()
  
  # Plot PCA-based UMAP
  return(seurat.subset)
}

# Usage example
seurat.subset <- seurat_subclustering(seurat_obj, 
                                      n_variable_features = 2000, 
                                      pca_dims = 20, 
                                      ident_column = "seurat_clusters", 
                                      idents_to_keep = c(0,1,2,3), 
                                      invert = FALSE)
