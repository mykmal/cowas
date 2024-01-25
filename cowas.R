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
    help = "Path to a genotype matrix of cis-variants for the first protein. This should be
                a tab-separated file with a header line followed by one line per individual,
                containing individual IDs in the first column and variants coded as 0/1/2 in
                the remaining columns. [required]"
  ),
  make_option(
    c("--genotypes_b"),
    help = "Path to a genotype matrix of cis-variants for the second protein, formatted the
                same as genotypes_a. [required]"
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
    c("--r2_threshold"),
    default = 0.01,
    type = "double",
    help = "R^2 threshold for expression and co-expression prediction models. The COWAS
                association test will only be performed if all three models have a
                cross-validation R^2 value above this threshold. [default %default]"
  ),
  make_option(
    c("--cv_folds"),
    default = 10,
    type = "integer",
    help = "Number of cross-validation folds to use for training elastic net models and
                assessing their performance. [default %default]"
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

# Create a data table with all cis-variants for both proteins
genotypes <- cbind(genotypes_a, genotypes_b[, !"IID"])
duplicated_variants <- which(duplicated(names(genotypes)))
suppressWarnings(genotypes[, (duplicated_variants) := NULL])

# Free up memory, since we won't need to distinguish between variants for one protein vs another
rm(genotypes_a)
rm(genotypes_b)

# Remove allele codes after each rsid, in case files were created with `plink2 --recode A`
setnames(genotypes, gsub(pattern = "_.*", replacement = "", x = names(genotypes)))

# Fill in missing calls with the mode for each variant
# This is a reasonable imputation method when the missingness rate is very low
for (variant in names(genotypes)[-1]) {
  mode <- names(which.max(table(genotypes[[variant]])))
  set(x = genotypes, i = which(is.na(genotypes[[variant]])), j = variant, value = mode)
}

# If covariates were not provided, simply normalize the expression data
# Otherwise, also adjust each protein's expression for the provided covariates
if (is.na(opt$covariates)) {
  # Remove rows that contain NAs
  expression <-	na.omit(expression)
  
  # Subset to the common set of individuals
  individuals <- intersect(genotypes$IID, expression$IID)
  genotypes <- genotypes[IID %in% individuals, ]
  expression <- expression[IID %in% individuals, ]
  
  # Normalize expression levels for each protein
  for (protein in c(opt$protein_a, opt$protein_b)) {
    expression_vector <- expression[[protein]]
    expression[, (protein) := scale(expression_vector)]
  }
} else {
  covariates <- fread(file = opt$covariates, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
  
  # Remove rows that contain NAs
  expression <- na.omit(expression)
  covariates <- na.omit(covariates)
  
  # Subset to the common set of individuals (including covariates in the intersection)
  individuals <- Reduce(intersect,
                        list(genotypes$IID,
                             expression$IID,
                             covariates$IID))
  genotypes <- genotypes[IID %in% individuals, ]
  expression <- expression[IID %in% individuals, ]
  covariates <- covariates[IID %in% individuals, ]
  
  # Adjust for covariates and normalize
  # After this, we check if expression levels are constant or have NAs
  for (protein in c(opt$protein_a, opt$protein_b)) {
    regression <- summary(lm(expression[[protein]] ~ ., data = covariates[, !"IID"]))
    expression[, (protein) := scale(regression$residuals)]
    
    if (anyNA(expression[[protein]]) || var(expression[[protein]]) <= 0) {
      stop("Expression of ", protein, " is either constant or NA. Skipping the pair ", opt$protein_a, " - ", opt$protein_b, ".")
    }
  }
}

# Load variant data and harmonize it with genotype data -------------------------------------------

# Read variants and remove duplicates from overlapping cis-regions
snps <- fread(file = opt$snps, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
              select = c("ID", "REF", "ALT"))
snps <- snps[ID != "ID", ]
snps <- snps[!duplicated(ID), ]

# Subset genotypes to those variants present in the annotation matrix
keep <- c(1, which(names(genotypes) %in% snps$ID))
genotypes <- genotypes[, ..keep]

# Normalize the genotype data for each variant, and then remove the SNP if it's monomorphic
for (rsid in names(genotypes)[-1]) {
  set(x = genotypes, j = rsid, value = scale(genotypes[[rsid]]))
  if (anyNA(genotypes[[rsid]]) || var(genotypes[[rsid]]) <= 0) {
    genotypes[, (rsid) := NULL]
  }
}

# Update the variant annotation matrix
snps <- snps[ID %in% names(genotypes), ]

# Load GWAS data and flip alleles -----------------------------------------------------------------

# Load GWAS summary statistics and merge them with the variant annotation matrix
gwas <- fread(file = opt$gwas, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
              select = c("variant_id", "effect_allele", "other_allele", "z_score", "n_cases", "n_controls"))
gwas <- na.omit(gwas)
gwas <- merge(gwas, snps, by.x = "variant_id", by.y = "ID", sort = FALSE)

# Flip the GWAS summary statistics when necessary
gwas <- gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]
gwas[effect_allele != ALT, flipped := TRUE]
gwas[flipped == TRUE, `:=` (effect_allele_flipped = other_allele,
                            other_allele_flipped = effect_allele,
                            z_score = -1 * as.numeric(z_score))]
gwas[flipped == TRUE, `:=` (effect_allele = effect_allele_flipped,
                            other_allele = other_allele_flipped)]
gwas[, c("REF", "ALT", "flipped", "effect_allele_flipped", "other_allele_flipped") := NULL]

# Remove any genotyped variants not present in the gwas
keep <- c(1, which(names(genotypes) %in% gwas$variant_id))
genotypes <- genotypes[, ..keep]

# Train and evaluate the prediction models --------------------------------------------------------

# Match up the genotype and expression data by sample ID, since we need to remove the IID colummn before model training
data_merged <- merge(expression, genotypes, by = "IID")
protein_columns <- c(opt$protein_a, opt$protein_b)
expression <- data_merged[, ..protein_columns]
rsids <- names(genotypes)[-1]
genotypes <- data_merged[, ..rsids]

# Free up memory
rm(data_merged)

# The glmnet package only accepts matrices and vectors as input
expression <- as.matrix(expression)
genotypes <- as.matrix(genotypes)

# Fit an elastic net model for protein_a
model_a <- cv.glmnet(x = genotypes, y = expression[, 1],
                     family = "gaussian", type.measure = "mse",
                     alpha = 0.5, nfolds = opt$cv_folds,
                     standardize = FALSE, intercept = TRUE)

# Fit an elastic net model for protein_b
model_b <- cv.glmnet(x = genotypes, y = expression[, 2],
                     family = "gaussian", type.measure = "mse",
                     alpha = 0.5, nfolds = opt$cv_folds,
                     standardize = FALSE, intercept = TRUE)

# Calculate the correlation between the expression of protein_a and protein_b,
# conditional on cis-variants for the genes coding both proteins
expression_a_imputed <- predict(model_a, newx = genotypes, s = "lambda.1se", type = "response")
expression_b_imputed <- predict(model_b, newx = genotypes, s = "lambda.1se", type = "response")
coexpression <- scale((expression[, 1] - expression_a_imputed) * (expression[, 2] - expression_b_imputed))[, 1]

# Fit an elastic net model for the co-expression of protein_a and protein_b
model_co <- cv.glmnet(x = genotypes, y = coexpression,
                      family = "gaussian", type.measure = "mse",
                      alpha = 0.5, nfolds = opt$cv_folds,
                      standardize = FALSE, intercept = TRUE)

# Save model weights
# The intercept is ignored because it will be numerically zero
model_a_weights <- coef(model_a, s = "lambda.1se")[-1, ]
model_b_weights <- coef(model_b, s = "lambda.1se")[-1, ]
model_co_weights <- coef(model_co, s = "lambda.1se")[-1, ]

# Check that all three models have at least some nonzero weights
if (var(model_a_weights) <= 0 || var(model_b_weights) <= 0 || var(model_co_weights) <= 0) {
  stop("At least one of the models in the pair ", opt$protein_a, " - ", opt$protein_b, " has no nonzero weights. This pair will be skipped.")
}

# Check that all three models meet the R^2 threshold
# R^2 = 1 - MSE / Var(y) but here Var(y) = 1
r2_a <- 1 - model_a$cvm[which(model_a$lambda == model_a$lambda.1se)]
r2_b <- 1 - model_b$cvm[which(model_b$lambda == model_b$lambda.1se)]
r2_co <- 1 - model_co$cvm[which(model_co$lambda == model_co$lambda.1se)]
if (r2_a < opt$r2_threshold || r2_b < opt$r2_threshold || r2_co < opt$r2_threshold) {
  stop("At least one of the models in the pair ", opt$protein_a, " - ", opt$protein_b, " fails the R^2 threshold. This pair will be skipped.")
}

# Test for association between imputed co-expression and disease ----------------------------------

# Get the sample sizes we are working with
n_genotypes <- nrow(genotypes)
n_gwas <- median(gwas$n_cases) + median(gwas$n_controls)

# Compute the correlation between variants and disease
genotype_disease_correlation <- gwas$z_score / sqrt(n_gwas - 2 + gwas$z_score^2)

# Formulas for computing the stage 2 effect size and its variance
stage2 <- function(qtl_weights) {
  # These formulas are derived from ordinary least squares (OLS) regression
  ld_matrix <- t(genotypes) %*% genotypes / n_genotypes
  product <- t(qtl_weights) %*% ld_matrix %*% qtl_weights
  
  # Check if the product matrix is invertible
  if (abs(det(product)) < 1e-20) {
    return(list(theta = 0,
                variance_theta = 0))
  }
  
  product_inverted <- solve(product)
  
  # See our paper for the derivation of the following three expressions -- it's too long to explain here
  theta <- product_inverted %*% t(qtl_weights) %*% genotype_disease_correlation
   
  stage2_residual_variance <- (1 / (n_gwas - ncol(qtl_weights))) * 
                                  (n_gwas - 1 - 2 * n_gwas * t(genotype_disease_correlation) %*% qtl_weights %*% theta +
                                      n_gwas * t(theta) %*% product %*% theta)
  
  variance_theta <- product_inverted %*% stage2_residual_variance / n_gwas
  
  return(list(theta = theta,
              variance_theta = variance_theta))
}

# Compute effect sizes and variances for each model
stage2_a <- stage2(qtl_weights = as.matrix(model_a_weights))
stage2_b <- stage2(qtl_weights = as.matrix(model_b_weights))
stage2_co <- stage2(qtl_weights = as.matrix(model_co_weights))
stage2_abc <- stage2(qtl_weights = cbind(model_a_weights, model_b_weights, model_co_weights))

# The "direct" variables refer to stage 2 models with only one term
if (stage2_a$variance_theta != 0) {
  wald_direct_a <- stage2_a$theta^2 / stage2_a$variance_theta
  pvalue_direct_a <- pchisq(wald_direct_a, 1, lower.tail = FALSE)
} else {
  pvalue_direct_a <- NA
}

if (stage2_b$variance_theta != 0) {
  wald_direct_b <- stage2_b$theta^2 / stage2_b$variance_theta
  pvalue_direct_b <- pchisq(wald_direct_b, 1, lower.tail = FALSE)
} else {
  pvalue_direct_b <- NA
}

if (stage2_co$variance_theta != 0) {
  wald_direct_co <- stage2_co$theta^2 / stage2_co$variance_theta
  pvalue_direct_co <- pchisq(wald_direct_co, 1, lower.tail = FALSE)
} else {
  pvalue_direct_co <- NA
}

# The "full" variables refer to the stage 2 model with three terms
if (all(stage2_abc$variance_theta != 0)) {
  wald_full_a <- stage2_abc$theta[1]^2 / stage2_abc$variance_theta[1]
  pvalue_full_a <- pchisq(wald_full_a, 1, lower.tail = FALSE)
  
  wald_full_b <- stage2_abc$theta[2]^2 / stage2_abc$variance_theta[2]
  pvalue_full_b <- pchisq(wald_full_b, 1, lower.tail = FALSE)
  
  wald_full_co <- stage2_abc$theta[3]^2 / stage2_abc$variance_theta[3]
  pvalue_full_co <- pchisq(wald_full_co, 1, lower.tail = FALSE)
} else {
  pvalue_full_a <- pvalue_full_b <- pvalue_full_co <- NA
}

# Append results to the output file
output <- c(opt$protein_a, opt$protein_b, n_genotypes, n_gwas,
            sum(model_a_weights != 0), r2_a,
            sum(model_b_weights != 0), r2_b,
            sum(model_co_weights != 0), r2_co,
            stage2_a$theta, stage2_a$variance_theta, pvalue_direct_a,
            stage2_b$theta, stage2_b$variance_theta, pvalue_direct_b,
            stage2_co$theta, stage2_co$variance_theta, pvalue_direct_co,
            stage2_abc$theta[1], stage2_abc$variance_theta[1], pvalue_full_a,
            stage2_abc$theta[2], stage2_abc$variance_theta[2], pvalue_full_b,
            stage2_abc$theta[3], stage2_abc$variance_theta[3], pvalue_full_co)

write.table(output, file = opt$out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)

