# Load necessary library
library(msigdbr)
library(Seurat)

# Function to get GOBP gene sets with at least 10 genes in the query
get_filtered_genesets <- function(msigdbr_list, query_genes, prefix="GOBP", min_genes=10) {
  # Filter gene sets by prefix (Select GOBP or HALLMARK or otherelse)
  filtered_genesets <- msigdbr_list[grep(prefix, names(msigdbr_list))]
  
  # Keep only gene sets with at least min_genes present in query_genes
  valid_genesets <- sapply(filtered_genesets, function(geneset) {
    length(intersect(geneset, query_genes)) >= min_genes
  })
  
  return(names(filtered_genesets[valid_genesets]))
}

# Get MSigDB human gene sets
msigdbr_human <- msigdbr(species = "Homo sapiens")
msigdbr_list <- split(msigdbr_human$gene_symbol, msigdbr_human$gs_name)

# Query genes
query_genes <- rownames(seurat_obj)

# Get filtered GOBP gene sets
GOBP_genesets <- get_filtered_genesets(msigdbr_list, query_genes)

# Add module scores to the Seurat object
seurat_obj <- AddModuleScore(
  object = seurat_obj, 
  features = msigdbr_list[GOBP_genesets], 
  name = "GOBP_pathway"
)
