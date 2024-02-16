#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(glmnet))

# Create command-line options ---------------------------------------------------------------------

option_list <- list(
  make_option(
    c("--protein_a"),
    help = "Name (or identifier) of the first protein in the co-expression pair. [required]"
  ),
  make_option(
    c("--protein_b"),
    help = "Name (or identifier) of the second protein in the co-expression pair. [required]"
  ),
  make_option(
    c("--genotypes_a"),
    help = "Path to a genotype matrix of variants to use as predictors for the first protein.
	            This should be a tab-separated file with a header line followed by one line
				per individual, containing individual IDs in the first column and variants
				with ALT alleles coded as 0..2 in the remaining columns. [required]"
  ),
  make_option(
    c("--genotypes_b"),
    help = "Path to a genotype matrix of variants to use as predictors for the second protein,
	            formatted the same as genotypes_a. [required]"
  ),
  make_option(
    c("--snps"),
    help = "Path to a PLINK 2.0 pvar file containing information on all of the variants
                present in --genotypes_a and --genotypes_b. [required]"
  ),
  make_option(
    c("--expression"),
    help = "Path to a matrix of expression levels for both proteins in the co-expression pair.
                This should be a tab-separated file with a header line followed by one line
                per individual, containing individual IDs in the first column and protein
                expression values in the remaining columns. [required]"
  ),
  make_option(
    c("--covariates"),
    help = "Path to a matrix of expression covariates. If specified, this should be a
                tab-separated file with a header line followed by one line per individual,
                containing individual IDs in the first column and covariates in the remaining
                columns. Note that categorical variables need to already be coded as dummy
                variables. [optional]"
  ),
    make_option(
    c("--gwas"),
    help = "Path to a GWAS summary statistics file for the disease of interest. This should be
                a tab-separated file with a header line followed by one line per variant,
                containing the following columns: variant_id, effect_allele, other_allele,
                z_score, n_cases, n_controls. Other columns are allowed but will be ignored.
                [required]"
  ),
  make_option(
    c("--out"),
    default = "output/ukb_proteins_all.tsv",
    help = "Path to a file where COWAS results will be stored. If the specified file already
                exists, a new line of results will be appended to its end. [default `%default`]"
  ),
  make_option(
    c("--model"),
	default = "stepwise",
	help = "The type of model to fit. Valid options are `stepwise` (linear regression with
	            both-direction stepwise variable selection), `ridge` (linear regression with
				a ridge penalty), and `elastic_net` (linear regression with an elastic net
				penalty). [default `%default`]"
  ),
  make_option(
    c("--r2_threshold"),
    default = 0.001,
    type = "double",
    help = "R^2 threshold for expression and co-expression prediction models. The COWAS
                association test will only be performed if all three models have a
                cross-validation R^2 value above this threshold. [default %default]"
  ),
  make_option(
    c("--cv_folds"),
    default = 10,
    type = "integer",
    help = "Number of cross-validation folds to use for training expression imputation models and
                assessing their predictive performance. Note that for ridge and elastic net models,
				the same folds are used to both tune the lambda hyperparameter and estimate R^2.
				Stepwise regression has no hyperparameters to tune, so cross-validation is only used
				to estimate R^2. [default %default]"
  ),
  make_option(
    c("--rank_normalize"),
	action = "store_true",
	default = TRUE,
	help = "Perform inverse-rank normalization (aka quantile normalization) on the expression
	            phenotypes before fitting models. If FALSE, expression values will simply be
				centered and scaled. [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (anyNA(opt[-7])) {
  stop("Some required parameters are missing. Run `cowas.R --help` for usage info.")
}

# Load genotype and expression data, normalize, and adjust for covariates -------------------------

genotypes_a <- fread(file = opt$genotypes_a, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
genotypes_b <- fread(file = opt$genotypes_b, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
expression <- fread(file = opt$expression, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
                    select = c("IID", opt$protein_a, opt$protein_b))

# Remove allele codes after each rsid, in case files were created with `plink2 --recode A`
setnames(genotypes_a, gsub(pattern = "_.*", replacement = "", x = names(genotypes_a)))
setnames(genotypes_b, gsub(pattern = "_.*", replacement = "", x = names(genotypes_b)))

# Create a data table with all predictor variants for both proteins
genotypes <- cbind(genotypes_a, genotypes_b[, !"IID"])
duplicated_variants <- which(duplicated(names(genotypes)))
suppressWarnings(genotypes[, (duplicated_variants) := NULL])

# Free up memory
variants_a <- names(genotypes_a)[-1]
variants_b <- names(genotypes_b)[-1]
rm(genotypes_a, genotypes_b)

# Subset to the common set of individuals with no missing values
expression <- na.omit(expression)
if (is.na(opt$covariates)) {
  individuals <- intersect(genotypes$IID, expression$IID)
  genotypes <- genotypes[IID %in% individuals, ]
  expression <- expression[IID %in% individuals, ]
} else {
  covariates <- fread(file = opt$covariates, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
  
  # Center and scale each covariate
  for (covariate in names(covariates)[-1]) {
    covariate_vector <- covariates[[covariate]]
    covariates[, (covariate) := scale(covariate_vector)]
  }
  covariates <- na.omit(covariates)
  
  individuals <- Reduce(intersect,
                        list(genotypes$IID,
                             expression$IID,
                             covariates$IID))
  genotypes <- genotypes[IID %in% individuals, ]
  expression <- expression[IID %in% individuals, ]
  covariates <- covariates[IID %in% individuals, ]
}

# Save the expression reference panel sample size
n_expression <- nrow(expression)

# Process each of the two proteins
for (protein in c(opt$protein_a, opt$protein_b)) {
  expression_vector <- expression[[protein]]
  
  # Perform inverse-rank normalization if requested
  if (opt$rank_normalize == TRUE) {
    # This offset corresponds to the commonly-used Blom transform
	offset <- 0.375
	
	# Compute the rank of each observation
	ranks <- rank(expression_vector, ties.method = "average")
	
	# Perform an inverse-rank normal transformation
	expression[, (protein) := stats::qnorm((ranks - offset) / (n_expression - 2 * offset + 1))]
  } else {
    # Otherwise, simply center and scale
	expression[, (protein) := scale(expression_vector)]
  }
  
  # Free up memory
  rm(expression_vector)
  
  # If covariates were provided, regress them out
  if (!is.na(opt$covariates)) {
    regression <- stats::lm(expression[[protein]] ~ ., data = covariates[, !"IID"])
	expression[, (protein) := regression$residuals]
  }
  
  # Finally, check if the processed expression levels are constant or have NAs
  if (anyNA(expression[[protein]]) || var(expression[[protein]]) <= 0) {
    stop("Expression of ", protein, " is either constant or NA after normalization and/or covariate adjustment. Skipping the pair ", opt$protein_a, " - ", opt$protein_b, ".")
  }
}

# Load variant data and harmonize it with genotype data -------------------------------------------

# Read variant list and remove duplicates
snps <- fread(file = opt$snps, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
              select = c("ID", "REF", "ALT"))
snps <- snps[ID != "ID", ]
snps <- snps[!duplicated(ID), ]

# Subset genotypes to those variants present in the variant list
keep <- c(1, which(names(genotypes) %in% snps$ID))
genotypes <- genotypes[, ..keep]

# Fill in missing calls with the mode for each variant
# This is a reasonable imputation method when the missingness rate is very low
for (variant in names(genotypes)[-1]) {
  mode <- names(which.max(table(genotypes[[variant]])))
  set(x = genotypes, i = which(is.na(genotypes[[variant]])), j = variant, value = mode)
}

# Normalize the genotype data for each variant, and then remove the SNP if it's monomorphic
for (rsid in names(genotypes)[-1]) {
  set(x = genotypes, j = rsid, value = scale(genotypes[[rsid]]))
  if (anyNA(genotypes[[rsid]]) || var(genotypes[[rsid]]) <= 0) {
    genotypes[, (rsid) := NULL]
  }
}

# Update the variant list
snps <- snps[ID %in% names(genotypes), ]

# Load GWAS data and flip alleles -----------------------------------------------------------------

# Load GWAS summary statistics and merge them with the variant list
gwas <- fread(file = opt$gwas, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
              select = c("variant_id", "effect_allele", "other_allele", "z_score", "n_cases", "n_controls"))
gwas <- na.omit(gwas)
gwas <- merge(gwas, snps, by.x = "variant_id", by.y = "ID", sort = FALSE)

# Flip the GWAS summary statistics when necessary
gwas <- gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]
gwas[effect_allele != ALT, flip := TRUE]
gwas[flip == TRUE, `:=` (effect_allele_flipped = other_allele,
                         other_allele_flipped = effect_allele,
                         z_score = -1 * as.numeric(z_score))]
gwas[flip == TRUE, `:=` (effect_allele = effect_allele_flipped,
                         other_allele = other_allele_flipped)]
gwas[, c("REF", "ALT", "flip", "effect_allele_flipped", "other_allele_flipped") := NULL]

# Remove any genotyped variants not present in the gwas
keep <- c(1, which(names(genotypes) %in% gwas$variant_id))
genotypes <- genotypes[, ..keep]

# Make sure at least two variants remain for each protein
variants_a <- variants_a[which(variants_a %in% names(genotypes))]
variants_b <- variants_b[which(variants_b %in% names(genotypes))]
if (length(variants_a) < 2 || length(variants_b) < 2) {
  stop("Fewer than two variants available for prediction. Skipping the pair ", opt$protein_a, " - ", opt$protein_b, ".")
}

# Define prediction models ------------------------------------------------------------------------

TrainStepwise <- function(z_a, z_b, z_both, x) {
  
  # Shuffle the data to break correlations before defining folds
  shuffled <- sample(n_expression)
  z_a <- z_a[shuffled, ]
  z_b <- z_b[shuffled, ]
  z_both <- z_both[shuffled, ]
  x <- x[shuffled, ]
  
  # lm() is easier to use with outcome and predictor variables in one data table
  data_a <- cbind(x[, 1], z_a)
  setnames(data_a, 1, "p_a")
  data_b <- cbind(x[, 2], z_b)
  setnames(data_b, 1, "p_b")
  
  # Free up memory
  rm(shuffled, z_a, z_b, x)
  
  # Define folds for cross-validation
  folds <- cut(seq(1, n_expression), breaks = opt$cv_folds, labels = FALSE)
  
  # Data table for storing imputed expression and co-expression in left-out folds
  imputed <- setDT(as.data.frame(matrix(data = NA, nrow = n_expression, ncol = 3)))
  setnames(imputed, c("a", "b", "co"))
  
  # The full model formulas
  formula_full_a <- as.formula(paste0("p_a ~ ", paste0(names(data_a)[-1], collapse = " + ")))
  formula_full_b <- as.formula(paste0("p_b ~ ", paste0(names(data_b)[-1], collapse = " + ")))
  formula_full_co <- as.formula(paste0("x_co ~ ", paste0(names(z_both), collapse = " + ")))
  
  for (k in 1:opt$cv_folds) {
    test_indices <- which(folds == k, arr.ind = TRUE)
	
	# Training data (k-1 folds)
	data_a_train <- data_a[!test_indices, ]
	data_b_train <- data_b[!test_indices, ]
	z_both_train <- z_both[!test_indices, ]
	
	# Test data (1 fold)
	data_a_test <- data_a[test_indices, ]
	data_b_test <- data_b[test_indices, ]
	z_both_test <- z_both[test_indices, ]
	
	# The starting models
	lm_null_a <- stats::lm(p_a ~ 1, data = data_a_train)
	lm_null_b <- stats::lm(p_b ~ 1, data = data_b_train)
	
	# Both-direction stepwise selection with AIC as the model selection criterion
	step_a <- stats::step(lm_null_a, scope = formula_full_a, direction = "both", trace = 0)
	step_b <- stats::step(lm_null_b, scope = formula_full_b, direction = "both", trace = 0)
	
	# Compute conditional co-expression
	x_a_imputed <- fitted(step_a)
	x_b_imputed <- fitted(step_b)
	x_co <- scale((data_a_train$p_a - x_a_imputed) * (data_b_train$p_b - x_b_imputed))[, 1]
	
	# Fit a model to predict co-expression
	data_co_train <- cbind(x_co, z_both_train)
	lm_null_co <- stats::lm(x_co ~ 1, data = data_co_train)
	step_co <- stats::step(lm_null_co, scope = formula_full_co, direction = "both", trace = 0)
	
	# Make predictions in the left-out fold
	imputed[test_indices, "a"] <- predict(step_a, data_a_test, type = "response")
	imputed[test_indices, "b"] <- predict(step_b, data_b_test, type = "response")
	imputed[test_indices, "co"] <- predict(step_co, z_both_test, type = "response")
  }
  
  # Fit full-sample models for each protein
  lm_null_a <- stats::lm(p_a ~ 1, data = data_a)
  lm_null_b <- stats::lm(p_b ~ 1, data = data_b)
  step_a <- stats::step(lm_null_a, scope = formula_full_a, direction = "both", trace = 0)
  step_b <- stats::step(lm_null_b, scope = formula_full_b, direction = "both", trace = 0)
  
  # Compute full-sample co-expression
  x_a_imputed <- fitted(step_a)
  x_b_imputed <- fitted(step_b)
  x_co <- scale((data_a$p_a - x_a_imputed) * (data_b$p_b - x_b_imputed))[, 1]
  
  # Fit a full-sample model for co-expression
  data_co <- cbind(x_co, z_both)
  lm_null_co <- stats::lm(x_co ~ 1, data = data_co)
  step_co <- stats::step(lm_null_co, scope = formula_full_co, direction = "both", trace = 0)
  
  # Save model weights
  # The intercept is ignored because it will be numerically zero
  model_a_weights <- coef(step_a)[-1]
  model_b_weights <- coef(step_b)[-1]
  model_co_weights <- coef(step_co)[-1]
  
  # Compute R^2 on left-out folds
  # R^2 = 1 - MSE / Var(y) = 1 - MSE
  # because here Var(y) = 1
  r2_a <- 1 - mean((data_a$p_a - imputed$a)^2)
  r2_b <- 1 - mean((data_b$p_b - imputed$b)^2)
  r2_co <- 1 - mean((fitted(step_co) - imputed$co)^2)
  
  return(list(weights_a = model_a_weights, weights_b = model_b_weights, weights_co = model_co_weights,
              r2_a = r2_a, r2_b = r2_b, r2_co = r2_co))
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
                       alpha = alpha, nfolds = opt$cv_folds,
                       standardize = FALSE, intercept = TRUE)
  
  # Fit an elastic net model for protein_b
  model_b <- cv.glmnet(x = z_b, y = x[, 2],
                       family = "gaussian", type.measure = "mse",
                       alpha = alpha, nfolds = opt$cv_folds,
                       standardize = FALSE, intercept = TRUE)
  
  # Calculate the correlation between the expression of protein_a and protein_b,
  # conditional on predictor variants for both proteins
  x_a_imputed <- predict(model_a, newx = z_a, s = "lambda.min", type = "response")
  x_b_imputed <- predict(model_b, newx = z_b, s = "lambda.min", type = "response")
  x_co <- scale((x[, 1] - x_a_imputed) * (x[, 2] - x_b_imputed))[, 1]
  
  # Fit an elastic net model for the co-expression of protein_a and protein_b
  model_co <- cv.glmnet(x = z_both, y = x_co,
                        family = "gaussian", type.measure = "mse",
                        alpha = alpha, nfolds = opt$cv_folds,
                        standardize = FALSE, intercept = TRUE)
  
  # Save model weights
  # The intercept is ignored because it will be numerically zero
  model_a_weights <- coef(model_a, s = "lambda.min")[-1, ]
  model_b_weights <- coef(model_b, s = "lambda.min")[-1, ]
  model_co_weights <- coef(model_co, s = "lambda.min")[-1, ]
  
  # R^2 = 1 - MSE / Var(y) = 1 - MSE
  # because here Var(y) = 1
  r2_a <- 1 - model_a$cvm[which(model_a$lambda == model_a$lambda.min)]
  r2_b <- 1 - model_b$cvm[which(model_b$lambda == model_b$lambda.min)]
  r2_co <- 1 - model_co$cvm[which(model_co$lambda == model_co$lambda.min)]
  
  return(list(weights_a = model_a_weights, weights_b = model_b_weights, weights_co = model_co_weights,
              r2_a = r2_a, r2_b = r2_b, r2_co = r2_co))
}

# Train and evaluate the prediction models --------------------------------------------------------

# Match up the genotype and expression data by sample ID, since we need to remove the IID colummn before model training
data_merged <- merge(expression, genotypes, by = "IID")
protein_columns <- c(opt$protein_a, opt$protein_b)
expression <- data_merged[, ..protein_columns]
rsids <- names(genotypes)[-1]
genotypes <- data_merged[, ..rsids]

# Free up memory
rm(data_merged)

# Extract protein-specific variants from the joint genotype matrix
genotypes_a <- genotypes[, ..variants_a]
genotypes_b <- genotypes[, ..variants_b]

# Train the requested model type
if (opt$model == "stepwise") {
  trained_output <- TrainStepwise(genotypes_a, genotypes_b, genotypes, expression)
} else if (opt$model == "ridge") {
  trained_output <- TrainGlmnet(genotypes_a, genotypes_b, genotypes, expression, alpha = 0)
} else if (opt$model == "elastic_net") {
  trained_output <- TrainGlmnet(genotypes_a, genotypes_b, genotypes, expression, alpha = 0.5)
}

# Check that all three models have at least some variation in weights among variants
if (var(trained_output$weights_a) <= 0 || var(trained_output$weights_b) <= 0 || var(trained_output$weights_co) <= 0) {
  stop("At least one of the models in the pair ", opt$protein_a, " - ", opt$protein_b, " has no nonzero weights. This pair will be skipped.")
}

# Check that all three models meet the R^2 threshold
if (trained_output$r2_a < opt$r2_threshold || trained_output$r2_b < opt$r2_threshold || trained_output$r2_co < opt$r2_threshold) {
  stop("At least one of the models in the pair ", opt$protein_a, " - ", opt$protein_b, " fails the R^2 threshold. This pair will be skipped.")
}

# Test for association between imputed co-expression and disease ----------------------------------

# Approximate the GWAS sample size
n_gwas <- median(gwas$n_cases) + median(gwas$n_controls)

# Compute the correlation between variants and disease
genotype_disease_correlation <- gwas$z_score / sqrt(n_gwas - 2 + gwas$z_score^2)

# Convert from data table to matrix so that we can perform matrix operations
genotypes <- as.matrix(genotypes)
genotype_disease_correlation <- as.matrix(genotype_disease_correlation)

# Compute the LD matrix
ld_matrix <- t(genotypes) %*% genotypes / n_expression

# Formulas for computing the stage 2 effect size and its variance
stage2 <- function(qtl_weights) {
  # This formula is derived from ordinary least squares (OLS) regression
  product <- t(qtl_weights) %*% ld_matrix %*% qtl_weights
  
  # Check if the product matrix is invertible
  if (abs(det(product)) < 1e-20) {
    return(list(theta = NA,
                variance_theta = NA,
				rss = NA))
  }
  
  product_inverted <- solve(product)
  
  # See our paper for the derivation of the following three expressions -- it's too long to explain here
  theta <- product_inverted %*% t(qtl_weights) %*% genotype_disease_correlation
  
  rss <- n_gwas * (1 - 2 * t(genotype_disease_correlation) %*% qtl_weights %*% theta + t(theta) %*% product %*% theta) - 1
  rss <- as.numeric(rss)
  
  variance_theta <- product_inverted * rss / ((n_gwas - ncol(qtl_weights) - 1) * n_gwas)
  
  return(list(theta = theta,
              variance_theta = variance_theta,
			  rss = rss))
}

# Fill in zeros for missing variants in the single-protein models to ensure that dimensions match the co-expression model
weights_padded_a <- weights_padded_b <- numeric(length(trained_output$weights_co))
names(weights_padded_a) <- names(weights_padded_b) <- names(trained_output$weights_co)
weights_padded_a[names(trained_output$weights_a)] <- trained_output$weights_a
weights_padded_b[names(trained_output$weights_b)] <- trained_output$weights_b

# Compute effect sizes and variances for each model
stage2_a <- stage2(as.matrix(weights_padded_a))
stage2_b <- stage2(as.matrix(weights_padded_b))
stage2_co <- stage2(as.matrix(trained_output$weights_co))
stage2_abc <- stage2(cbind(weights_padded_a, weights_padded_b, trained_output$weights_co))

# The "direct" variables refer to stage 2 models with only one term
if (!is.na(stage2_a$variance_theta)) {
  wald_direct_a <- as.numeric(stage2_a$theta)^2 / as.numeric(stage2_a$variance_theta)
  pvalue_direct_a <- pchisq(wald_direct_a, 1, lower.tail = FALSE)
} else {
  pvalue_direct_a <- NA
}

if (!is.na(stage2_b$variance_theta)) {
  wald_direct_b <- as.numeric(stage2_b$theta)^2 / as.numeric(stage2_b$variance_theta)
  pvalue_direct_b <- pchisq(wald_direct_b, 1, lower.tail = FALSE)
} else {
  pvalue_direct_b <- NA
}

if (!is.na(stage2_co$variance_theta)) {
  wald_direct_co <- as.numeric(stage2_co$theta)^2 / as.numeric(stage2_co$variance_theta)
  pvalue_direct_co <- pchisq(wald_direct_co, 1, lower.tail = FALSE)
} else {
  pvalue_direct_co <- NA
}

# The "full" variables refer to the stage 2 model with three terms
if (!anyNA(stage2_abc$variance_theta)) {
  wald_full_a <- stage2_abc$theta[1,1]^2 / stage2_abc$variance_theta[1,1]
  pvalue_full_a <- pchisq(wald_full_a, 1, lower.tail = FALSE)
  
  wald_full_b <- stage2_abc$theta[2,1]^2 / stage2_abc$variance_theta[2,2]
  pvalue_full_b <- pchisq(wald_full_b, 1, lower.tail = FALSE)
  
  wald_full_co <- stage2_abc$theta[3,1]^2 / stage2_abc$variance_theta[3,3]
  pvalue_full_co <- pchisq(wald_full_co, 1, lower.tail = FALSE)
} else {
  pvalue_full_a <- pvalue_full_b <- pvalue_full_co <- NA
}

# Perform an F-test to determine if the stage 2 model with three terms is significantly better than a null model
f_statistic <- ((n_gwas - 4) * (n_gwas - 1 - as.numeric(stage2_abc$rss))) / (3 * as.numeric(stage2_abc$rss))
f_pvalue <- pf(f_statistic, 3, n_gwas - 4, lower.tail = FALSE)

# Append results to the output file
output <- c(opt$protein_a, opt$protein_b, n_expression, n_gwas,
            sum(trained_output$weights_a != 0), trained_output$r2_a,
            sum(trained_output$weights_b != 0), trained_output$r2_b,
            sum(trained_output$weights_co != 0), trained_output$r2_co,
            stage2_a$theta, stage2_a$variance_theta, pvalue_direct_a,
            stage2_b$theta, stage2_b$variance_theta, pvalue_direct_b,
            stage2_co$theta, stage2_co$variance_theta, pvalue_direct_co,
            stage2_abc$theta[1], stage2_abc$variance_theta[1], pvalue_full_a,
            stage2_abc$theta[2], stage2_abc$variance_theta[2], pvalue_full_b,
            stage2_abc$theta[3], stage2_abc$variance_theta[3], pvalue_full_co,
			f_statistic, f_pvalue)

write.table(t(output), file = opt$out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)

