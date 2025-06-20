#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))

# Create command-line options ---------------------------------------------------------------------

option_list <- list(
  make_option(
    c("--protein_a"),
    action = "store",
    type = "character",
    default = NA,
    help = "Name (or identifier) of the first protein in the co-expression pair. [required]"
  ),
  make_option(
    c("--protein_b"),
    action = "store",
    type = "character",
    default = NA,
    help = "Name (or identifier) of the second protein in the co-expression pair. [required]"
  ),
  make_option(
    c("--genotypes_a"),
    action = "store",
    type = "character",
    default = NA,
    help = "Path to a genotype matrix of variants to use as predictors for the first protein.
            This should be a tab-separated file with a header line followed by one line per
            individual, containing individual IDs in the first column and variants with effect
            allele dosages coded as 0..2 in the remaining columns. [required]"
  ),
  make_option(
    c("--genotypes_b"),
    action = "store",
    type = "character",
    default = NA,
    help = "Path to a genotype matrix of variants to use as predictors for the second protein,
            formatted the same as --genotypes_a. [required]"
  ),
  make_option(
    c("--expression"),
    action = "store",
    type = "character",
    default = NA,
    help = "Path to a matrix of expression levels for both proteins in the co-expression pair.
            This should be a tab-separated file with a header line followed by one line per
            individual, containing individual IDs in the first column and protein expression
            values in the remaining columns. [required]"
  ),
  make_option(
    c("--covariates"),
    action = "store",
    type = "character",
    default = NA,
    help = "Path to a matrix of expression covariates. If specified, this should be a
            tab-separated file with a header line followed by one line per individual,
            containing individual IDs in the first column and covariates in the remaining
            columns. Note that categorical variables need to already be coded as dummy
            variables. [optional]"
  ),
  make_option(
    c("--out_folder"),
    action = "store",
    type = "character",
    default = "cowas_weights",
    help = "Path to a folder where COWAS weights will be stored. A TSV file with performance
            metrics for each model will also be saved here. [default: '%default']"
  ),
  make_option(
    c("--model"),
    action = "store",
    type = "character",
    default = "elastic_net",
    help = "The type of model to fit. Valid options are 'stepwise' (linear regression with
            both-direction stepwise variable selection by AIC), 'ridge' (linear regression
            with an L2 penalty), 'lasso' (linear regression with an L1 penalty), and
            'elastic_net' (linear regression with a linear combination of the L1 and L2
            penalties). [default: '%default']"
  ),
  make_option(
    c("--cores"),
    action = "store",
    type = "integer",
    default = 1L,
    help = "Number of cores to use for parallelization. The default value disables parallel
            computation. [default: %default]"
  ),
  make_option(
    c("--cor_threshold"),
    action = "store",
    type = "double",
    default = 0.03,
    help = "Correlation threshold for expression and co-expression imputation models. COWAS
            model weights will only be saved if the Pearson correlation between measured
            and predicted expression (calculated on a held-out test set) is above this
            threshold for all three models. [default: %default]"
  ),
  make_option(
    c("--rank_normalize"),
    action = "store",
    type = "logical",
    default = TRUE,
    help = "Perform a rank-based inverse-normal transformation on the expression phenotypes
            before fitting models. If FALSE, expression values will simply be centered
            and scaled. [default: %default]"
  ),
  make_option(
    c("--product_based"),
    action = "store",
    type = "logical",
    default = FALSE,
    help = "Setting this parameter to TRUE will train models for product-based COWAS. In the
            default version of COWAS, the outcome of the co-expression model is the product
            of single-protein model residuals. In the product-based version of COWAS, the
            outcome of the co-expression model is instead the product of observed expression
            levels. See our paper for details. [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (anyNA(opt[-6])) {
  stop("Some required parameters are missing. Run './cowas_train.R --help' for usage info.")
}

# Load the glmnet package if penalized regression is requested
if (opt$model == "ridge" || opt$model == "lasso" || opt$model == "elastic_net") {
  suppressMessages(library(glmnet))
}

# Enable parallel computation if requested
if (opt$cores > 1) {
  suppressMessages(library(doMC))
  registerDoMC(cores = opt$cores)
  
  setDTthreads(threads = opt$cores, restore_after_fork = TRUE)
  use_cores <- TRUE
} else {
  use_cores <- FALSE
}

# Load genotype, expression, and covariate data ---------------------------------------------------

genotypes_a <- fread(file = opt$genotypes_a, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
genotypes_b <- fread(file = opt$genotypes_b, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
expression <- fread(file = opt$expression, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
                    select = c("IID", opt$protein_a, opt$protein_b))

# Remove allele codes after each rsid, in case files were created with 'plink2 --recode A'
setnames(genotypes_a, gsub(pattern = "_.*", replacement = "", x = names(genotypes_a)))
setnames(genotypes_b, gsub(pattern = "_.*", replacement = "", x = names(genotypes_b)))

# Create a data table with all predictor variants for both proteins
variants_a <- names(genotypes_a)[-1]
variants_b <- names(genotypes_b)[-1]
genotypes <- genotypes_a[genotypes_b[, -which(names(genotypes_b) %in% variants_a), with = FALSE], on = "IID"]
rm(genotypes_a, genotypes_b)

# Subset to the common set of individuals with no missing values
expression <- na.omit(expression)
if (is.na(opt$covariates)) {
  individuals <- intersect(genotypes$IID, expression$IID)
  genotypes <- genotypes[IID %in% individuals, ]
  expression <- expression[IID %in% individuals, ]
  
  # Sort each dataset by IID so that they match
  setorder(genotypes, IID)
  setorder(expression, IID)
  
} else {
  covariates <- fread(file = opt$covariates, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
  covariates <- na.omit(covariates)
  
  individuals <- Reduce(intersect,
                        list(genotypes$IID,
                             expression$IID,
                             covariates$IID))
  genotypes <- genotypes[IID %in% individuals, ]
  expression <- expression[IID %in% individuals, ]
  covariates <- covariates[IID %in% individuals, ]
  
  # Sort each dataset by IID so that they match
  setorder(genotypes, IID)
  setorder(expression, IID)
  setorder(covariates, IID)
  
  # Center and scale each covariate
  for (column in names(covariates)[-1]) {
    set(x = covariates, j = column, value = scale(covariates[[column]]))
    
    # Remove the covariate if it has NAs or is constant for this subset of individuals
    if (anyNA(covariates[[column]]) || var(covariates[[column]]) <= 0) {
      covariates[, (column) := NULL]
    }
  }
}

# Save the expression reference panel sample size
n_expression <- nrow(expression)

# Normalize, adjust expression for covariates, and impute missing genotype calls ------------------

# Process each of the two proteins
for (protein in c(opt$protein_a, opt$protein_b)) {
  
  # Perform rank-based inverse-normal transformation if requested
  if (opt$rank_normalize == TRUE) {
    # This offset corresponds to the commonly-used Blom transform
    offset <- 0.375
    
    # Compute the rank of each observation
    ranks <- rank(expression[[protein]], ties.method = "average")
    
    # Perform the transformation
    expression[, (protein) := stats::qnorm((ranks - offset) / (n_expression - 2 * offset + 1))]
  } else {
    # Otherwise, simply center and scale
    set(x = expression, j = protein, value = scale(expression[[protein]]))
  }
  
  # If covariates were provided, regress them out
  if (!is.na(opt$covariates)) {
    regression <- stats::lm(expression[[protein]] ~ ., data = covariates[, !"IID"])
    expression[, (protein) := regression$residuals]
  }
  
  # Check if the processed expression levels are constant or have NAs
  if (anyNA(expression[[protein]]) || var(expression[[protein]]) <= 0) {
    stop("Expression of ", protein, " is either constant or NA after normalization and/or covariate adjustment. Skipping the pair ", opt$protein_a, "_", opt$protein_b, ".")
  }
}

# Free up memory
rm(covariates)

# Fill in missing genotype calls with the mode for each variant.
# This is a reasonable imputation method when the missingness rate is very low.
for (rsid in names(genotypes)[-1]) {
  
  mode <- names(which.max(table(genotypes[[rsid]])))
  set(x = genotypes, i = which(is.na(genotypes[[rsid]])), j = rsid, value = mode)
  
  # Standardize the genotype data for this variant, and then remove it if it's monomorphic
  set(x = genotypes, j = rsid, value = scale(genotypes[[rsid]]))
  if (anyNA(genotypes[[rsid]]) || var(genotypes[[rsid]]) <= 0) {
    genotypes[, (rsid) := NULL]
  }
}

# Make sure at least two variants remain for each protein
variants_a <- variants_a[which(variants_a %in% names(genotypes))]
variants_b <- variants_b[which(variants_b %in% names(genotypes))]
if (length(variants_a) < 2 || length(variants_b) < 2) {
  stop("Fewer than two variants available for prediction. Skipping the pair ", opt$protein_a, "_", opt$protein_b, ".")
}

# Define prediction models ------------------------------------------------------------------------

TrainStepwise <- function(z_a, z_b, z_both, x) {
  
  # lm() is easier to use with outcome and predictor variables in one data table
  data_a <- cbind(x[, 1], z_a)
  setnames(data_a, 1, "x_a")
  data_b <- cbind(x[, 2], z_b)
  setnames(data_b, 1, "x_b")
  
  # Free up memory
  rm(z_a, z_b, x)
  
  # The full model formulas
  formula_full_a <- as.formula(paste0("x_a ~ ", paste0(names(data_a)[-1], collapse = " + ")))
  formula_full_b <- as.formula(paste0("x_b ~ ", paste0(names(data_b)[-1], collapse = " + ")))
  
  # Fit models for each protein
  lm_null_a <- stats::lm(x_a ~ 1, data = data_a)
  lm_null_b <- stats::lm(x_b ~ 1, data = data_b)
  step_a <- stats::step(lm_null_a, scope = formula_full_a, direction = "both", trace = 0)
  step_b <- stats::step(lm_null_b, scope = formula_full_b, direction = "both", trace = 0)
  
  # Estimate the expected conditional correlation if running the standard version of COWAS,
  # otherwise compute the product of observed expression levels.
  if (opt$product_based == FALSE) {
    pred_a <- predict(object = step_a, newdata = data_a, type = "response")
    pred_b <- predict(object = step_b, newdata = data_b, type = "response")
    x_co <- (data_a$x_a - pred_a) * (data_b$x_b - pred_b)
    rm(pred_a, pred_b)
  } else {
    x_co <- data_a$x_a * data_b$x_b
  }
  
  # Create a data table containing co-expression values and all variants
  data_co <- cbind(x_co, z_both)
  rm(x_co, z_both)
  
  # Train a model to predict the two proteins' co-expression
  formula_full_co <- as.formula(paste0("x_co ~ ", paste0(names(data_co)[-1], collapse = " + ")))
  lm_null_co <- stats::lm(x_co ~ 1, data = data_co)
  step_co <- stats::step(lm_null_co, scope = formula_full_co, direction = "both", trace = 0)
  
  # Store fitted model weights.
  # The intercept is ignored because it only captures non-genetic effects.
  model_a_weights <- coef(step_a)[-1]
  model_b_weights <- coef(step_b)[-1]
  model_co_weights <- coef(step_co)[-1]
  
  return(list(weights_a = model_a_weights, weights_b = model_b_weights, weights_co = model_co_weights,
              model_a = step_a, model_b = step_b, model_co = step_co))
}

TrainGlmnet <- function(z_a, z_b, z_both, x, alpha) {
  
  # The glmnet package only accepts matrices and vectors as input
  z_a <- as.matrix(z_a)
  z_b <- as.matrix(z_b)
  z_both <- as.matrix(z_both)
  x <- as.matrix(x)
  
  # Fit an elastic net model for protein_a
  model_a <- cv.glmnet(x = z_a, y = x[, 1],
                       family = "gaussian", type.measure = "mse",
                       alpha = alpha, nfolds = 10,
                       standardize = FALSE, intercept = TRUE,
                       parallel = use_cores)
  
  # Fit an elastic net model for protein_b
  model_b <- cv.glmnet(x = z_b, y = x[, 2],
                       family = "gaussian", type.measure = "mse",
                       alpha = alpha, nfolds = 10,
                       standardize = FALSE, intercept = TRUE,
                       parallel = use_cores)
  
  # Estimate the expected conditional correlation if running the standard version of COWAS,
  # otherwise compute the product of observed expression levels.
  if (opt$product_based == FALSE) {
    pred_a <- predict(object = model_a, newx = z_a, s = "lambda.min", type = "response")
    pred_b <- predict(object = model_b, newx = z_b, s = "lambda.min", type = "response")
    coex <- (x[, 1] - pred_a) * (x[, 2] - pred_b)
    rm(pred_a, pred_b)
  } else {
    coex <- x[, 1] * x[, 2]
  }
  
  # Fit an elastic net model to predict the two proteins' co-expression
  model_co <- cv.glmnet(x = z_both, y = coex,
                        family = "gaussian", type.measure = "mse",
                        alpha = alpha, nfolds = 10,
                        standardize = FALSE, intercept = TRUE,
                        parallel = use_cores)
  
  # Store fitted model weights.
  # The intercept is ignored because it only captures non-genetic effects.
  model_a_weights <- coef(model_a, s = "lambda.min")[-1, ]
  model_b_weights <- coef(model_b, s = "lambda.min")[-1, ]
  model_co_weights <- coef(model_co, s = "lambda.min")[-1, ]
  
  return(list(weights_a = model_a_weights, weights_b = model_b_weights, weights_co = model_co_weights,
              model_a = model_a, model_b = model_b, model_co = model_co))
}

# Train and evaluate the prediction models --------------------------------------------------------

# Make sure the genotype and expression data are matched by sample ID, since we need to remove the IID column
expression <- expression[match(genotypes$IID, IID), ]
expression[, IID := NULL]
genotypes[, IID := NULL]

# Extract protein-specific variants from the joint genotype matrix
genotypes_a <- genotypes[, ..variants_a]
genotypes_b <- genotypes[, ..variants_b]

# Split data into training and test subsets
test_indices <- sample(x = n_expression, size = floor(0.2 * n_expression))

# Train the requested model type on the training set
if (opt$model == "stepwise") {
  training_output <- TrainStepwise(genotypes_a[!test_indices, ], genotypes_b[!test_indices, ], genotypes[!test_indices, ],
                                   expression[!test_indices, ])
} else if (opt$model == "ridge") {
  training_output <- TrainGlmnet(genotypes_a[!test_indices, ], genotypes_b[!test_indices, ], genotypes[!test_indices, ],
                                 expression[!test_indices, ], alpha = 0)
} else if (opt$model == "lasso") {
  training_output <- TrainGlmnet(genotypes_a[!test_indices, ], genotypes_b[!test_indices, ], genotypes[!test_indices, ],
                                 expression[!test_indices, ], alpha = 1)
} else if (opt$model == "elastic_net") {
  training_output <- TrainGlmnet(genotypes_a[!test_indices, ], genotypes_b[!test_indices, ], genotypes[!test_indices, ],
                                 expression[!test_indices, ], alpha = 0.5)
}

# Check that each of the three models has at least one variant with a nonzero weight
if (sum(training_output$weights_a != 0) <= 0 || sum(training_output$weights_b != 0) <= 0 || sum(training_output$weights_co != 0) <= 0) {
  stop("At least one of the models in the pair ", opt$protein_a, "_", opt$protein_b, " has no nonzero weights. This pair will be skipped.")
}

# Get predicted values on the test set
if (opt$model == "stepwise") {
  imputed_test_a <- predict(object = training_output$model_a, newdata = genotypes_a[test_indices, ],
                            type = "response")
  imputed_test_b <- predict(object = training_output$model_b, newdata = genotypes_b[test_indices, ],
                            type = "response")
  imputed_test_co <- predict(object = training_output$model_co, newdata = genotypes[test_indices, ],
                             type = "response")
} else {
  imputed_test_a <- predict(object = training_output$model_a, newx = as.matrix(genotypes_a[test_indices, ]),
                            s = "lambda.min", type = "response")
  imputed_test_b <- predict(object = training_output$model_b, newx = as.matrix(genotypes_b[test_indices, ]),
                            s = "lambda.min", type = "response")
  imputed_test_co <- predict(object = training_output$model_co, newx = as.matrix(genotypes[test_indices, ]),
                             s = "lambda.min", type = "response")
}

# Check that the predicted values are not constant
if (sd(imputed_test_a) <= 0 || sd(imputed_test_b) <= 0 || sd(imputed_test_co) <= 0) {
  stop("Imputed expression or co-expression is constant for at least one of the models in the pair ", opt$protein_a, "_", opt$protein_b, ". This pair will be skipped.")
}

# Fit full-sample models
rm(training_output)
if (opt$model == "stepwise") {
  full_output <- TrainStepwise(genotypes_a, genotypes_b, genotypes,
                               expression)
} else if (opt$model == "ridge") {
  full_output <- TrainGlmnet(genotypes_a, genotypes_b, genotypes,
                             expression, alpha = 0)
} else if (opt$model == "lasso") {
  full_output <- TrainGlmnet(genotypes_a, genotypes_b, genotypes,
                             expression, alpha = 1)
} else if (opt$model == "elastic_net") {
  full_output <- TrainGlmnet(genotypes_a, genotypes_b, genotypes,
                             expression, alpha = 0.5)
}

# Get the number of variants with nonzero weights in each model
n_nonzero_a <- sum(full_output$weights_a != 0)
n_nonzero_b <- sum(full_output$weights_b != 0)
n_nonzero_co <- sum(full_output$weights_co != 0)

# Check that each of the three models has at least one variant with a nonzero weight
if (n_nonzero_a <= 0 || n_nonzero_b <= 0 || n_nonzero_co <= 0) {
  stop("At least one of the models in the pair ", opt$protein_a, "_", opt$protein_b, " has no nonzero weights. This pair will be skipped.")
}

# Get predicted values on the full data set
if (opt$model == "stepwise") {
  imputed_full_a <- predict(object = full_output$model_a, newdata = genotypes_a,
                            type = "response")
  imputed_full_b <- predict(object = full_output$model_b, newdata = genotypes_b,
                            type = "response")
} else {
  imputed_full_a <- predict(object = full_output$model_a, newx = as.matrix(genotypes_a),
                            s = "lambda.min", type = "response")
  imputed_full_b <- predict(object = full_output$model_b, newx = as.matrix(genotypes_b),
                            s = "lambda.min", type = "response")
}

# If running the standard version of COWAS, estimate the conditional covariance for all samples
if (opt$product_based == FALSE) {
  set(x = expression, j = "coexpression", value = (expression[[opt$protein_a]] - imputed_full_a) * (expression[[opt$protein_b]] - imputed_full_b))
} else {
  # If instead running the product-based version of COWAS, calculate the product of the two proteins' observed expression levels
  set(x = expression, j = "coexpression", value = expression[[opt$protein_a]] * expression[[opt$protein_b]])
}

# Check that the predicted expression values and estimated co-expression values are not constant
if (sd(imputed_full_a) <= 0 || sd(imputed_full_b) <= 0 || sd(expression$coexpression) <= 0) {
  stop("Imputed expression or co-expression is constant for at least one of the models in the pair ", opt$protein_a, "_", opt$protein_b, ". This pair will be skipped.")
}

# For each model, compute the correlation and an association P value on the test set
cor_test_a <- stats::cor.test(expression[test_indices, ][[opt$protein_a]], imputed_test_a,
                              alternative = "two.sided", method = "pearson")
cor_test_b <- stats::cor.test(expression[test_indices, ][[opt$protein_b]], imputed_test_b,
                              alternative = "two.sided", method = "pearson")
cor_test_co <- stats::cor.test(expression[test_indices, ][["coexpression"]], imputed_test_co,
                               alternative = "two.sided", method = "pearson")

# Check that all three models pass the correlation threshold
if (cor_test_a$estimate < opt$cor_threshold || cor_test_b$estimate < opt$cor_threshold || cor_test_co$estimate < opt$cor_threshold) {
  stop("At least one of the models in the pair ", opt$protein_a, "_", opt$protein_b, " fails the correlation threshold. This pair will be skipped.")
}

# Save model weights and performance metrics ------------------------------------------------------

# Weights for all three models are saved within a single list in an RDS file
weights_final <- full_output[c("weights_a", "weights_b", "weights_co")]
saveRDS(weights_final, file = paste0(opt$out_folder, "/", opt$protein_a, "_", opt$protein_b, ".weights.rds"))

# Gather performance metrics
metrics <- c(opt$protein_a, opt$protein_b, n_expression,
             n_nonzero_a, cor_test_a$estimate, cor_test_a$p.value, (cor_test_a$estimate)^2,
             n_nonzero_b, cor_test_b$estimate, cor_test_b$p.value, (cor_test_b$estimate)^2,
             n_nonzero_co, cor_test_co$estimate, cor_test_co$p.value, (cor_test_co$estimate)^2)

# Append the performance metrics to the output file
write.table(t(metrics), file = paste0(opt$out_folder, "/performance_metrics.tsv"),
            quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE, append = TRUE)

