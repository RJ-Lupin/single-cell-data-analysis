# Reference of Source code: "https://github.com/immunogenomics/masc"
# Citation: "Fonseka, Chamith Y., et al. "Mixed-effects association of single cells identifies an expanded effector CD4+ T cell subset in rheumatoid arthritis." Science translational medicine 10.463 (2018): eaaq0305."

# The following code are copied from "https://github.com/immunogenomics/masc/blob/master/R/masc.R" and "https://github.com/immunogenomics/masc/blob/master/R/utils.R"

#' MASC - Mixed effect modeling of Associations of Single Cells
#'
#' @param dataset A data frame containing the contrast factor, random, and fixed effects for the model
#' @param cluster A factor indicating cluster assignments for each cell
#' @param contrast A vector indicating the variable to be tested for association with cluster abundance. Must match a column in dataset.
#' @param random_effects A vector indicating which terms should be modeled as random effects covariates. Terms listed must match columns in dataset.
#' @param fixed_effects A vector indicating which terms should be modeled as fixed effects covariates. Terms listed must match columns in dataset.
#' @param save_models Should MASC save the mixed-effects model objects generated for each cluster?
#' @param save_model_dir Location to save mixed-effect model objects. Defaults to current working directory.
#' @param verbose TRUE/FALSE
#'
#' @return data frame containing calculated association p-values and odds ratios for each cluster tested

MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                 verbose = FALSE, save_models = FALSE, save_model_dir = NULL) {
  # Check inputs
  if (is.factor(dataset[[contrast]]) == FALSE) {
    stop("Specified contrast term is not coded as a factor in dataset")
  }

  # Generate design matrix from cluster assignments
  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)

  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]

  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + ")),
                        collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(fixed_effects, collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      # For now, do not allow models without mixed effects terms
      stop("No random effects specified")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
    model_rhs <- "1" # only includes intercept
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }

  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]

  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                   model_rhs), collapse = ""))
    # Run null and full mixed-effects models
    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci <- confint.merMod(full_model, method = "Wald",
                                  parm = contrast_lvl2)
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }

  # Organize results into output dataframe
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                       size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))

  # Return MASC results and save models if specified
  if (save_models == TRUE) {
    saveModelObj(cluster_models, save_dir = save_model_dir)
    return(output)
  } else {
    return(output)
  }
}

saveModelObj <- function(model_obj, save_name = "masc.modelobj.rds", save_dir = NULL) {
  # If directory is unspecified, use current working directory
  if(is.null(save_dir)) {
    message("save_model_dir is unspecified, saving model objects to current directory")
    save_dir <- getwd()
  }
  # Try to save file unless it already exists
  if(file.exists(save_name) == FALSE) {
    saveRDS(model_obj, file = save_name)
    message(paste("Models saved to", file.path(save_dir, save_name)))
  } else {
    warning(paste(save_name, "already exists in directory, did not overwrite"))
  }
}

# Usage
library(lme4)

masc.input = seurat.obj@meta.data
masc.input$cell_type = factor(masc.input$cell_type)
masc.input$contrast_group = factor(masc.input$contrast_group)
masc.input$sample_id = factor(masc.input$sample_id)

# run masc
masc_res = MASC(dataset = masc.input,
                cluster = masc.input$cell_type,
                contrast = "contrast_group",
                random_effects = "sample_id",
                verbose = TRUE, 
                save_models = F)

masc_res$FDR = p.adjust(masc_res$model.pvalue, method = "fdr", n = length(masc_res$model.pvalue))
masc_res$cell_type = substring(masc_res$cluster,8,)
masc_res$log2OR = log2(masc_res[,grep("OR$",colnames(masc_res))])
masc_res$log2ci.lower = log2(masc_res[,grep("ci.lower",colnames(masc_res))])
masc_res$log2ci.upper = log2(masc_res[,grep("ci.upper",colnames(masc_res))])
masc_res$log10_pval = (-log10(masc_res$model.pvalue))
masc_res$log10_p.adj = (-log10(masc_res$FDR))
masc_res$cell_type = factor(masc_res$cell_type, levels=rev(levels(masc.input$cell_type)))

# Visualization
ggplot(masc_res, aes(x=log2OR, y=cell_type)) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = 'dashed') +
  geom_errorbarh(aes(xmax = log2ci.upper, xmin = log2ci.lower), size = .5, height = .2, color = 'gray50') +
  geom_point(aes(size = log10_p.adj))

