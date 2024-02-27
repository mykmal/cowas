#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))

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
				with effect allele dosages coded as 0..2 in the remaining columns. [required]"
  ),
  make_option(
    c("--genotypes_b"),
    help = "Path to a genotype matrix of variants to use as predictors for the second protein,
	            formatted the same as --genotypes_a. [required]"
  ),
  make_option(
    c("--snps"),
    help = "Path to a table specifying the reference and effect alleles for each variant present
				in --genotypes_a and --genotypes_b. This should be a tab-separated file with a
				header line followed by one line per variant, containing the following three
				columns: ID, REF, ALT. Other columns are allowed but will be ignored. [required]"
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
                n_cases, n_controls, z_score. Other columns are allowed but will be ignored.
                [required]"
  ),
  make_option(
    c("--out"),
    default = "output/results.tsv",
    help = "Path to a file where COWAS results will be stored. If the specified file already
                exists, a new line of results will be appended to its end. [default `%default`]"
  ),
  make_option(
    c("--model"),
	default = "stepwise",
	help = "The type of model to fit. Valid options are `stepwise` (linear regression with
	            both-direction stepwise variable selection), `ridge` (linear regression with
				a ridge penalty), `lasso` (linear regression with a lasso penalty), and
				`elastic_net` (linear regression with an elastic net penalty). [default `%default`]"
  ),
  make_option(
    c("--cores"),
	default = "1",
	type = "integer",
	help = "Number of cores to use for parallelization. The default value disables parallel
				computation. [default `%default`]"
  ),
  make_option(
    c("--r2_threshold"),
    default = 0.001,
    type = "double",
    help = "R^2 threshold for expression and co-expression prediction models. The COWAS
                association test will only be performed if all three models have a predictive
                R^2 value (calculated on a held-out test set) above this threshold.
				[default %default]"
  ),
  make_option(
    c("--rank_normalize"),
	action = "store_true",
	default = TRUE,
	help = "Perform a rank-based inverse-normal transformation (aka quantile normalization)
				on the expression phenotypes before fitting models. If FALSE, expression values
				will simply be centered and scaled. [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (anyNA(opt[-7])) {
  stop("Some required parameters are missing. Run `cowas.R --help` for usage info.")
}

# Load the glmnet package if penalized regression is requested
if (opt$model == "ridge" || opt$model == "lasso" || opt$model == "elastic_net") {
  suppressMessages(library(glmnet))
}

# Enable parallel computation if requested
if (opt$cores > 1L) {
  suppressMessages(library(doMC))
  registerDoMC(cores = opt$cores)
  
  suppressMessages(setDTthreads(threads = opt$cores, restore_after_fork = TRUE))
  use_cores <- TRUE
} else {
  use_cores <- FALSE
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
  covariates <- na.omit(covariates)
  
  individuals <- Reduce(intersect,
                        list(genotypes$IID,
                             expression$IID,
                             covariates$IID))
  genotypes <- genotypes[IID %in% individuals, ]
  expression <- expression[IID %in% individuals, ]
  covariates <- covariates[IID %in% individuals, ]
  
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

# Process each of the two proteins
for (protein in c(opt$protein_a, opt$protein_b)) {
  
  # Perform quantile normalization if requested
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
    stop("Expression of ", protein, " is either constant or NA after normalization and/or covariate adjustment. Skipping the pair ", opt$protein_a, " - ", opt$protein_b, ".")
  }
}

# Finally, compute the co-expression between the two proteins
set(x = expression, j = "coexpression", value = expression[[opt$protein_a]] * expression[[opt$protein_b]])

# Load variant data and harmonize it with genotype data -------------------------------------------

# Read variant list and remove duplicates
snps <- fread(file = opt$snps, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
              select = c("ID", "REF", "ALT"))
snps <- snps[!duplicated(ID), ]

# Subset genotypes to those variants present in the variant list
keep <- c(1, which(names(genotypes) %in% snps$ID))
genotypes <- genotypes[, ..keep]

# Fill in missing calls with the mode for each variant
# This is a reasonable imputation method when the missingness rate is very low
for (rsid in names(genotypes)[-1]) {
  
  mode <- names(which.max(table(genotypes[[rsid]])))
  set(x = genotypes, i = which(is.na(genotypes[[rsid]])), j = rsid, value = mode)
  
  # Standardize the genotype data for this variant, and then remove it if it's monomorphic
  set(x = genotypes, j = rsid, value = scale(genotypes[[rsid]]))
  if (anyNA(genotypes[[rsid]]) || var(genotypes[[rsid]]) <= 0) {
    genotypes[, (rsid) := NULL]
  }
}

# Update the variant list
snps <- snps[ID %in% names(genotypes)[-1], ]

# Load GWAS data and flip alleles -----------------------------------------------------------------

# Load GWAS summary statistics and merge them with the variant list
gwas <- fread(file = opt$gwas, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
              select = c("variant_id", "effect_allele", "other_allele", "n_cases", "n_controls", "z_score"))
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
  
  # lm() is easier to use with outcome and predictor variables in one data table
  data_a <- cbind(x[, 1], z_a)
  setnames(data_a, 1, "x_a")
  data_b <- cbind(x[, 2], z_b)
  setnames(data_b, 1, "x_b")
  data_co <- cbind(x[, 3], z_both)
  setnames(data_co, 1, "x_co")
  
  # Free up memory
  rm(z_a, z_b, z_both, x)
  
  # The full model formulas
  formula_full_a <- as.formula(paste0("x_a ~ ", paste0(names(data_a)[-1], collapse = " + ")))
  formula_full_b <- as.formula(paste0("x_b ~ ", paste0(names(data_b)[-1], collapse = " + ")))
  formula_full_co <- as.formula(paste0("x_co ~ ", paste0(names(data_co)[-1], collapse = " + ")))
  
  # Fit models for each protein
  lm_null_a <- stats::lm(x_a ~ 1, data = data_a)
  lm_null_b <- stats::lm(x_b ~ 1, data = data_b)
  lm_null_co <- stats::lm(x_co ~ 1, data = data_co)
  step_a <- stats::step(lm_null_a, scope = formula_full_a, direction = "both", trace = 0)
  step_b <- stats::step(lm_null_b, scope = formula_full_b, direction = "both", trace = 0)
  step_co <- stats::step(lm_null_co, scope = formula_full_co, direction = "both", trace = 0)
  
  # Save fitted model weights
  # The intercept is ignored because it will be numerically zero
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
  
  # Fit an elastic net model for the co-expression of protein_a and protein_b
  model_co <- cv.glmnet(x = z_both, y = x[, 3],
                        family = "gaussian", type.measure = "mse",
                        alpha = alpha, nfolds = 10,
                        standardize = FALSE, intercept = TRUE,
						parallel = use_cores)
  
  # Save model weights
  # The intercept is ignored because it will be numerically zero
  model_a_weights <- coef(model_a, s = "lambda.1se")[-1, ]
  model_b_weights <- coef(model_b, s = "lambda.1se")[-1, ]
  model_co_weights <- coef(model_co, s = "lambda.1se")[-1, ]
  
  return(list(weights_a = model_a_weights, weights_b = model_b_weights, weights_co = model_co_weights,
              model_a = model_a, model_b = model_b, model_co = model_co))
}

# Train and evaluate the prediction models --------------------------------------------------------

# Match up the genotype and expression data by sample ID, since we need to remove the IID column before model training
data_merged <- merge(expression, genotypes, by = "IID")
expression_columns <- c(opt$protein_a, opt$protein_b, "coexpression")
expression <- data_merged[, ..expression_columns]
rsids <- names(genotypes)[-1]
genotypes <- data_merged[, ..rsids]

# Free up memory
rm(data_merged)

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

# Check that all three models have at least some variation in weights among variants
if (var(training_output$weights_a) <= 0 || var(training_output$weights_b) <= 0 || var(training_output$weights_co) <= 0) {
  stop("At least one of the models in the pair ", opt$protein_a, " - ", opt$protein_b, " has no nonzero weights. This pair will be skipped.")
}

# Get predicted values on the test set
imputed_a <- predict(training_output$model_a, genotypes_a[test_indices, ], type = "response")
imputed_b <- predict(training_output$model_b, genotypes_b[test_indices, ], type = "response")
imputed_co <- predict(training_output$model_co, genotypes[test_indices, ], type = "response")

# Compute R^2 on the test set
# R^2 = 1 - MSE / Var(y) = 1 - MSE
# because here Var(y) = 1
r2_a <- 1 - mean((expression[test_indices, 1] - imputed_a)^2)
r2_b <- 1 - mean((expression[test_indices, 2] - imputed_b)^2)
r2_co <- 1 - mean((expression[test_indices, 3] - imputed_co)^2)

# Check that all three models pass the R^2 threshold
if (r2_a < opt$r2_threshold || r2_b < opt$r2_threshold || r2_co < opt$r2_threshold) {
  stop("At least one of the models in the pair ", opt$protein_a, " - ", opt$protein_b, " fails the R^2 threshold. This pair will be skipped.")
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
weights_padded_a <- weights_padded_b <- numeric(length(full_output$weights_co))
names(weights_padded_a) <- names(weights_padded_b) <- names(full_output$weights_co)
weights_padded_a[names(full_output$weights_a)] <- full_output$weights_a
weights_padded_b[names(full_output$weights_b)] <- full_output$weights_b

# Compute effect sizes and variances for each model
stage2_a <- stage2(as.matrix(weights_padded_a))
stage2_b <- stage2(as.matrix(weights_padded_b))
stage2_co <- stage2(as.matrix(full_output$weights_co))
stage2_abc <- stage2(cbind(weights_padded_a, weights_padded_b, full_output$weights_co))

# Convert from 1x1 matrix to numeric
stage2_a$theta <- as.numeric(stage2_a$theta)
stage2_a$variance_theta <- as.numeric(stage2_a$variance_theta)
stage2_b$theta <- as.numeric(stage2_b$theta)
stage2_b$variance_theta <- as.numeric(stage2_b$variance_theta)
stage2_co$theta <- as.numeric(stage2_co$theta)
stage2_co$variance_theta <- as.numeric(stage2_co$variance_theta)

# The "direct" variables refer to stage 2 models with only one term
if (!is.na(stage2_a$variance_theta)) {
  wald_direct_a <- stage2_a$theta^2 / stage2_a$variance_theta
  pvalue_direct_a <- stats::pchisq(wald_direct_a, 1, lower.tail = FALSE)
} else {
  pvalue_direct_a <- NA
}

if (!is.na(stage2_b$variance_theta)) {
  wald_direct_b <- stage2_b$theta^2 / stage2_b$variance_theta
  pvalue_direct_b <- stats::pchisq(wald_direct_b, 1, lower.tail = FALSE)
} else {
  pvalue_direct_b <- NA
}

if (!is.na(stage2_co$variance_theta)) {
  wald_direct_co <- stage2_co$theta^2 / stage2_co$variance_theta
  pvalue_direct_co <- stats::pchisq(wald_direct_co, 1, lower.tail = FALSE)
} else {
  pvalue_direct_co <- NA
}

# The "full" variables refer to the stage 2 model with three terms
if (!anyNA(stage2_abc$variance_theta)) {
  wald_full_a <- stage2_abc$theta[1,1]^2 / stage2_abc$variance_theta[1,1]
  pvalue_full_a <- stats::pchisq(wald_full_a, 1, lower.tail = FALSE)
  
  wald_full_b <- stage2_abc$theta[2,1]^2 / stage2_abc$variance_theta[2,2]
  pvalue_full_b <- stats::pchisq(wald_full_b, 1, lower.tail = FALSE)
  
  wald_full_co <- stage2_abc$theta[3,1]^2 / stage2_abc$variance_theta[3,3]
  pvalue_full_co <- stats::pchisq(wald_full_co, 1, lower.tail = FALSE)
} else {
  pvalue_full_a <- pvalue_full_b <- pvalue_full_co <- NA
}

# Perform an F-test to determine if the stage 2 model with three terms is significantly better than a null model
f_statistic <- ((n_gwas - 4) * (n_gwas - 1 - as.numeric(stage2_abc$rss))) / (3 * as.numeric(stage2_abc$rss))
f_pvalue <- stats::pf(f_statistic, 3, n_gwas - 4, lower.tail = FALSE)

# Append results to the output file
output <- c(opt$protein_a, opt$protein_b, n_expression, n_gwas,
            sum(full_output$weights_a != 0), r2_a,
            sum(full_output$weights_b != 0), r2_b,
            sum(full_output$weights_co != 0), r2_co,
            stage2_a$theta, stage2_a$variance_theta, pvalue_direct_a,
            stage2_b$theta, stage2_b$variance_theta, pvalue_direct_b,
            stage2_co$theta, stage2_co$variance_theta, pvalue_direct_co,
            stage2_abc$theta[1,1], stage2_abc$variance_theta[1,1], pvalue_full_a,
            stage2_abc$theta[2,1], stage2_abc$variance_theta[2,2], pvalue_full_b,
            stage2_abc$theta[3,1], stage2_abc$variance_theta[3,3], pvalue_full_co,
			f_statistic, f_pvalue)

write.table(t(output), file = opt$out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)

