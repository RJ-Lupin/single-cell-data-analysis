library(Seurat) # https://satijalab.org/seurat/ (doi: 10.1016/j.cell.2021.04.048)
library(harmony) # https://github.com/immunogenomics/harmony (doi: 10.1038/s41592-019-0619-0)
library(dplyr)

# Function for data integration using the Harmony package, followed by Seurat's standard scRNA-seq analysis workflow
perform_seurat_harmony_analysis <- function(seurat_object, 
                                            seed = 12345,
                                            n_variable_features = 2000,
                                            pca_dims = 20,
                                            harmony_dims = 10,
                                            harmony_group = "sample", # A column named 'sample' in the metadata of the Seurat object, which indicates batch or sample information
                                            harmony_max_iter = 30) {
  set.seed(seed)
  
  # Standard Seurat workflow
  seurat <- seurat_object %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = n_variable_features) %>%
    ScaleData() %>%
    RunPCA()
  
  # Elbow plot for PCA
  print(ElbowPlot(seurat, ndims = 50))

  # optional process for Harmony integration
  seurat <- seurat %>%
    RunUMAP(dims = 1:pca_dims) %>%
    FindNeighbors(dims = 1:pca_dims) %>%
    FindClusters()
  
  # Harmony integration
  harmony_seurat <- RunHarmony(seurat, group.by = harmony_group, max_iter = harmony_max_iter)
  
  # Elbow plot for Harmony
  print(ElbowPlot(harmony_seurat, reduction = "harmony", ndims = 50))
  
  # Downstream analysis with Harmony results
  harmony_seurat <- harmony_seurat %>%
    RunUMAP(reduction = "harmony", dims = 1:harmony_dims) %>%
    FindNeighbors(reduction = "harmony", dims = 1:harmony_dims) %>%
    FindClusters()
  
  # Plot Harmony-integrated UMAP
  return(list(standard = seurat, harmony = harmony_seurat))
}

# usage
harmony_results = perform_seurat_harmony_analysis(seurat_object, n_variable_features = 2000, pca_dims = 20, harmony_dims = 10, harmony_group = "batch_key")
seurat_harmony = harmony_results$harmony
