#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))

# Create command-line options ---------------------------------------------------------------------

option_list <- list(
  make_option(
    c("--protein_a"),
    action = "store",
    type = "character",
    help = "Name (or identifier) of the first protein in the co-expression pair. [required]"
  ),
  make_option(
    c("--protein_b"),
    action = "store",
    type = "character",
    help = "Name (or identifier) of the second protein in the co-expression pair. [required]"
  ),
  make_option(
    c("--gwas"),
    action = "store",
    type = "character",
    help = "Path to a GWAS summary statistics file for the outcome trait of interest.
            This should be a tab-separated file with a header line followed by one line
            per variant, containing the following columns: variant_id, effect_allele,
            other_allele, z_score, n_samples. Other columns are allowed but will be
            ignored. Note that variant IDs must be consistent with those used in COWAS
            weights and in the LD reference genotypes. (In our provided weights,
            variants are denoted by rsID.) [required]"
  ),
  make_option(
    c("--weights"),
    action = "store",
    type = "character",
    default = "cowas_weights",
    help = "Path to a folder containing trained COWAS weights in RDS format, as saved by
            cowas_train.R. [default: '%default']"
  ),
  make_option(
    c("--alleles"),
    action = "store",
    type = "character",
    default = "cowas_model_alleles.tsv",
    help = "Path to a table specifying the reference and effect alleles for each variant
            present in the trained COWAS models. This should be a tab-separated file with
            a header line followed by one line per variant, containing the following three
            columns: ID, REF, ALT. Other columns are allowed but will be ignored.
            [default: '%default']"
  ),
  make_option(
    c("--ld_reference"),
    action = "store",
    type = "character",
    help = "Path to an individual-level genotype file to use for computing a linkage
            disequilibrium (LD) reference panel. This should be a tab-separated file with
            a header line followed by one line per individual, containing individual IDs
            in the first column and variants with allele dosages coded as 0..2 in the
            remaining columns. Variant IDs in the header must be consistent with those
            used in COWAS weights and in the GWAS. (In our provided weights, variants are
            denoted by rsID.) [required]"
  ),
  make_option(
    c("--out"),
    action = "store",
    type = "character",
    default = "cowas_results.tsv",
    help = "Path to a file where COWAS results will be saved. If the specified file
            already exists, a new line of results will be appended to its end.
            [default: '%default']"
  ),
  make_option(
    c("--cores"),
    action = "store",
    type = "integer",
    default = 1,
    help = "Number of cores to use for parallelization. The default value disables
            multi-threaded computation. [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (anyNA(opt)) {
  stop("Some required parameters are missing. Run 'cowas.R --help' for usage info.")
}

# Set the requested number of cores
setDTthreads(threads = opt$cores, restore_after_fork = TRUE)

# Load GWAS data, model weights, allele file, and LD reference genotypes --------------------------

# Read the GWAS summary statistics
gwas <- fread(file = opt$gwas, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
              select = c("variant_id", "effect_allele", "other_allele", "z_score", "n_samples"))
gwas <- na.omit(gwas)

# Read the trained model weights
weights <- readRDS(paste0(opt$weights, "/", opt$protein_a, "_", opt$protein_b, ".weights.rds"))
weights_a <- weights$weights_a
weights_b <- weights$weights_b
weights_co <- weights$weights_co
rm(weights)

# Read the allele list
snps <- fread(file = opt$alleles, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
              select = c("ID", "REF", "ALT"))

# Read the LD reference genotypes
genotypes <- fread(file = opt$ld_reference, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)

# Remove allele codes after each rsid, in case files were created with 'plink2 --recode A'
setnames(genotypes, gsub(pattern = "_.*", replacement = "", x = names(genotypes)))

# Subset all data sources to a common set of variants ---------------------------------------------

# Variants included in at least one model
model_variants <- unique(c(names(weights_a),
                           names(weights_b),
                           names(weights_co)))

# Variants common to GWAS, models, allele list, and LD reference
common_variants <- Reduce(intersect,
                          list(gwas$variant_id,
                               model_variants,
                               snps$ID,
                               names(genotypes)))

# Subset the GWAS
gwas <- gwas[variant_id %in% common_variants, ]
setorder(gwas, variant_id)

# Subset the model weights
weights_a <- weights_a[intersect(common_variants, names(weights_a))]
weights_a <- weights_a[sort(names(weights_a))]
weights_b <- weights_b[intersect(common_variants, names(weights_b))]
weights_b <- weights_b[sort(names(weights_b))]
weights_co <- weights_co[intersect(common_variants, names(weights_co))]
weights_co <- weights_co[sort(names(weights_co))]

# Subset the allele list
snps <- snps[ID %in% common_variants, ]
setorder(snps, ID)

# Subset the LD reference genotype data
keep <- which(names(genotypes) %in% common_variants)
genotypes <- genotypes[, ..keep]
setcolorder(genotypes, sort(names(genotypes)))

# Flip GWAS alleles -------------------------------------------------------------------------------

gwas <- merge(gwas, snps, by.x = "variant_id", by.y = "ID", sort = FALSE)

# Remove variants with mismatching alleles
matches <- gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]$variant_id
matches <- sort(matches)
if (length(matches) < nrow(gwas)) {
  gwas <- gwas[variant_id %in% matches, ]
  setorder(gwas, variant_id)
  
  weights_a <- weights_a[intersect(matches, names(weights_a))]
  weights_a <- weights_a[sort(names(weights_a))]
  weights_b <- weights_b[intersect(matches, names(weights_b))]
  weights_b <- weights_b[sort(names(weights_b))]
  weights_co <- weights_co[intersect(matches, names(weights_co))]
  weights_co <- weights_co[sort(names(weights_co))]
  
  keep <- which(names(genotypes) %in% matches)
  genotypes <- genotypes[, ..keep]
  setcolorder(genotypes, sort(names(genotypes)))
}

# Flip the GWAS z-scores when necessary
gwas[effect_allele != ALT, z_score := -1 * as.numeric(z_score)]
gwas[, c("effect_allele", "other_allele", "REF", "ALT") := NULL]

# Compute an LD reference matrix from normalized genotypes ----------------------------------------

# Fill in missing calls with the mode for each variant.
# This is a reasonable imputation method when the missingness rate is very low.
for (rsid in names(genotypes)) {
  
  mode <- names(which.max(table(genotypes[[rsid]])))
  set(x = genotypes, i = which(is.na(genotypes[[rsid]])), j = rsid, value = mode)
  
  # Standardize the genotypes for this variant, and then remove it if it's monomorphic
  set(x = genotypes, j = rsid, value = scale(genotypes[[rsid]]))
  if (anyNA(genotypes[[rsid]]) || var(genotypes[[rsid]]) <= 0) {
    genotypes[, (rsid) := NULL]
    gwas <- gwas[variant_id != rsid, ]
    weights_a <- weights_a[names(weights_a) != rsid]
    weights_b <- weights_b[names(weights_b) != rsid]
    weights_co <- weights_co[names(weights_co) != rsid]
  }
}

# Make sure at least one variant remains in each model
if (length(weights_a) < 1 || length(weights_b) < 1 || length(weights_co) < 1) {
  stop("No variants remain in the intersection of model weights, GWAS effects, and the LD reference. Skipping the pair ", opt$protein_a, " - ", opt$protein_b, ".")
}

# Approximate the GWAS sample size
n_gwas <- median(gwas$n_samples)

# Save the LD reference sample size
n_reference <- nrow(genotypes)

# Compute the pseudocorrelation between each variant and the outcome trait
genotype_trait_correlations <- gwas$z_score / sqrt(n_gwas - 1 + gwas$z_score^2)

# Convert from data table to matrix so that we can perform matrix operations
genotypes <- as.matrix(genotypes)
genotype_trait_correlations <- as.matrix(genotype_trait_correlations)

# Compute the LD matrix
ld_matrix <- t(genotypes) %*% genotypes / n_reference

# Compute association between imputed (co-)expression and the trait -------------------------------

# Formulas for computing PWAS/COWAS effect sizes and their variances.
# See our paper for derivations.
compute_effect_size <- function(qtl_weights, n_terms) {

  product <- t(qtl_weights) %*% ld_matrix %*% qtl_weights
  
  # Check if the product matrix is invertible
  if (abs(det(product)) < 1e-20) {
    return(list(theta = matrix(data = NA, nrow = ncol(qtl_weights), ncol = 1),
                variance_theta = matrix(data = NA, nrow = ncol(qtl_weights), ncol = ncol(qtl_weights)),
                rss = NA))
  }
  
  product_inverted <- solve(product)
  
  theta <- product_inverted %*% t(qtl_weights) %*% genotype_trait_correlations
  
  rss <- n_gwas * (1 - 2 * t(genotype_trait_correlations) %*% qtl_weights %*% theta + t(theta) %*% product %*% theta) - 1
  rss <- as.numeric(rss)
  
  variance_theta <- product_inverted * rss / (n_gwas * (n_gwas - n_terms))
  
  return(list(theta = theta,
              variance_theta = variance_theta,
              rss = rss))
}

# Fill in zeros for variants without any weights so that dimensions match across all models
weights_padded_a <- weights_padded_b <- weights_padded_co <- numeric(ncol(genotypes))
names(weights_padded_a) <- names(weights_padded_b) <- names(weights_padded_co) <- colnames(genotypes)
weights_padded_a[names(weights_a)] <- weights_a
weights_padded_b[names(weights_b)] <- weights_b
weights_padded_co[names(weights_co)] <- weights_co

# Compute effect sizes and variances for each model
stage2_a <- compute_effect_size(as.matrix(weights_padded_a), 2)
stage2_b <- compute_effect_size(as.matrix(weights_padded_b), 2)
stage2_co <- compute_effect_size(as.matrix(weights_padded_co), 2)
stage2_abc <- compute_effect_size(cbind(weights_padded_a, weights_padded_b, weights_padded_co), 4)

# Perform association tests -----------------------------------------------------------------------

# Convert from 1x1 matrix to numeric
stage2_a$theta <- as.numeric(stage2_a$theta)
stage2_a$variance_theta <- as.numeric(stage2_a$variance_theta)
stage2_b$theta <- as.numeric(stage2_b$theta)
stage2_b$variance_theta <- as.numeric(stage2_b$variance_theta)
stage2_co$theta <- as.numeric(stage2_co$theta)
stage2_co$variance_theta <- as.numeric(stage2_co$variance_theta)

# The "direct" variables refer to expression-trait association models with only one regressor
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

# The "full" variables refer to the COWAS expression-trait association model, which has three regressors
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

# Perform an F-test to determine if the COWAS association model is significantly better than a null model
if (!is.na(stage2_abc$rss)) {
  stage2_abc$rss <- as.numeric(stage2_abc$rss)
  f_statistic <- ((n_gwas - 4) / (3 * stage2_abc$rss)) * (n_gwas - 1 - stage2_abc$rss)
  f_pvalue <- stats::pf(f_statistic, 3, n_gwas - 4, lower.tail = FALSE)
} else {
  f_statistic <- f_pvalue <- NA
}

# Save effect sizes and test results --------------------------------------------------------------

# Append results to the output file
output <- c(opt$protein_a, opt$protein_b, n_reference, n_gwas,
            stage2_a$theta, sqrt(stage2_a$variance_theta), pvalue_direct_a,
            stage2_b$theta, sqrt(stage2_b$variance_theta), pvalue_direct_b,
            stage2_abc$theta[1,1], sqrt(stage2_abc$variance_theta[1,1]), pvalue_full_a,
            stage2_abc$theta[2,1], sqrt(stage2_abc$variance_theta[2,2]), pvalue_full_b,
            stage2_abc$theta[3,1], sqrt(stage2_abc$variance_theta[3,3]), pvalue_full_co,
            f_statistic, f_pvalue)

write.table(t(output), file = opt$out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)

