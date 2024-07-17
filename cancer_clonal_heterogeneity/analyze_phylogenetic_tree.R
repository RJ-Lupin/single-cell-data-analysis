library(Seurat)
library(numbat)
library(ape)
library(phytools)

analyze_phylogenetic_tree <- function(out_dir) {
  # Initialize Numbat object
  numbat_obj <- Numbat$new(out_dir = out_dir)
  
  # Access the maximum likelihood tree
  tree <- numbat_obj$treeML
  
  # Calculate total phylogenetic events
  total_events <- nrow(tree$edge) - Ntip(tree) + 1
  
  # Calculate subclonal events
  # Assuming the root node is the most recent common ancestor (MRCA)
  mrca_node <- Ntip(tree) + 1
  subclonal_events <- nrow(tree$edge) - mrca_node + 1
  
  # Calculate subclonality
  subclonality <- subclonal_events / total_events
  
  # Analyze phylogenetic signal sidedness
  # Create a dataframe with sample and side information
  sample_data <- data.frame(
    sample = tree$tip.label,
    side = ifelse(grepl("right", tree$tip.label), 1,
                  ifelse(grepl("left", tree$tip.label), 0, 0.5))
  )
  
  # Calculate phylogenetic signal using Pagel's lambda
  phylosig_result <- phylosig(tree, sample_data$side, method = "lambda", test = TRUE)
  
  # Print results
  print(paste("Total phylogenetic events:", total_events))
  print(paste("Subclonal events:", subclonal_events))
  print(paste("Subclonality:", subclonality))
  print("Phylogenetic signal result:")
  print(phylosig_result)
}

# Example usage:
# analyze_phylogenetic_tree("/data/path")
