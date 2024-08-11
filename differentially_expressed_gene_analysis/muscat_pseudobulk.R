# Load required libraries
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(Seurat)
library(scater)
library(SingleCellExperiment)

# Set output directory
output_dir <- "path/to/output/directory"

# Preprocess SingleCellExperiment object
sce <- as.SingleCellExperiment(seurat)
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
sce <- prepSCE(sce, kid = "cluster_id", sid = "sample_id", gid = "group_id")

# Aggregate data
pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# Check available assays
print(assayNames(pb))

# Create design and contrast matrices
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("A-B", levels = mm)

# Define DE methods
muscat_de_methods <- c("edgeR", "DESeq2", "limma-trend", "limma-voom")

# Run DE analysis for each method
for (de_method in muscat_de_methods) {
  # Run DS analysis
  res <- pbDS(pb, design = mm, contrast = contrast, method = de_method)
  
  # Generate DEG output
  deg_output <- resDS(sce, res, bind = "row")
  
  # Write results to CSV
  filename <- paste0("muscat_", de_method, "_results.csv")
  write.csv(deg_output, file.path(output_dir, filename))
}
