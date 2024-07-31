# Define the function
perform_tukey_hsd <- function(data, response_var, main_effect, covariate) {
  
  # Create the formula
  # Creates a formula object using the response variable, main effect, and covariate. This formula is used to specify the model for the ANOVA analysis.
  formula <- as.formula(paste(response_var, "~", main_effect, "+", covariate))
  
  # Create the ANOVA model
  aov_model <- aov(formula, data = data)
  
  # Perform the Tukey HSD test
  # Applies the Tukey Honest Significant Difference (HSD) test to the ANOVA model to find significant differences between group means.
  res <- TukeyHSD(aov_model)
  
  # Return the result
  return(res)
}

# Create example data frame (assign your data frame here)
df_test <- seurat_obj@meta.data[, c("contrast_group", "batch", "signature_score")]

# Call the function
result <- perform_tukey_hsd(data = df_test, response_var = "signature_score", main_effect = "contrast_group", covariate = "batch")

# Print the result
print(result$contrast_group)
