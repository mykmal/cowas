#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(glmnet))

option_list <- list(
  make_option(
    c("-a", "--gene_a"),
    help = "Name (e.g. Ensembl ID) of the first gene in the co-expression pair. [required]"
  ),
  make_option(
    c("-b", "--gene_b"),
    help = "Name (e.g. Ensembl ID) of the second gene in the co-expression pair. [required]"
  ),
  make_option(
    c("-w", "--genotypes_a"),
    help = "A genotype matrix of cis-SNPs for the first gene. This should be a tab-separated
                file with a header line followed by one line per individual, containing
                individual IDs in the first column and variants coded as 0/1/2 in the
                remaining columns. [required]"
  ),
  make_option(
    c("-x", "--genotypes_b"),
    help = "A genotype matrix of cis-SNPs for the second gene, formatted the same as
                genotypes_a. [required]"
  ),
  make_option(
    c("-e", "--expression"),
    help = "A matrix of expression levels for both genes in the co-expression pair. This
                should be a tab-separated file with a header line followed by one line per
                individual, containing individual IDs in the first column and gene
                expression values in the remaining columns. [required]"
  ),
  make_option(
    c("-c", "--covariates"),
    help = "A matrix of covariates to residualize the expression values. If specified, this
                should be a tab-separated file with a header line followed by one line per
                individual, containing individual IDs in the first column and covariates in
                the remaining columns. [optional]"
  ),
  make_option(
    c("-s", "--snps"),
    help = "A table with SNP information. This should be a tab-separated file with one line
                per SNP and no header line, containing the following six columns:
                chromosome code, rsID, position in centimorgans (not used), base-pair
                coordinate, REF allele, ALT allele. [required]"
  ),
  make_option(
    c("-o", "--out"),
    default = "output",
    help = "Directory where COWAS weights and cross-validated performance metrics will be
                saved. [default `%default`]"
  ),
  make_option(
    c("-p", "--p_threshold"),
    default = 0.01,
    type = "double",
    help = "P-value threshold for nominal significance of expression and co-expression
                prediction models. If all models significantly explain some variability in
                their response variables at this threshold, COWAS weights will be computed
                and saved. [default %default]"
  ),
  make_option(
    c("-m", "--elnet_alpha"),
    default = 0.5,
    type = "double",
    help = "Mixing parameter (alpha) for elastic net regularization. [default %default]"
  ),
  make_option(
    c("-k", "--cv_folds"),
    default = 5,
    type = "integer",
    help = "Number of cross-validation folds to use for assessing model performance. Set
                to zero to skip cross-validation. [default %default]"
  ),
  make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Show verbose output. [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))


# Imputation models ------------------------------------------------------------------------------

# Train elastic net models to predict expression and co-expression
TrainCowasModels <- function(z_a, z_b, z_both, x_a, x_b) {
  
  z_a <- as.matrix(z_a[, !"IID"])
  z_b <- as.matrix(z_b[, !"IID"])
  z_both <- as.matrix(z_both[, !"IID"])
  
  # Fit an elastic net model for gene_a
  model_a <- cv.glmnet(x = z_a, y = x_a,
                       family = "gaussian", type.measure = "mse",
                       alpha = opt$elnet_alpha, nfolds = 5, intercept = FALSE)
  
  # Fit an elastic net model for gene_b
  model_b <- cv.glmnet(x = z_b, y = x_b,
                       family = "gaussian", type.measure = "mse",
                       alpha = opt$elnet_alpha, nfolds = 5, intercept = FALSE)
  
  # Calculate the correlation between the expression of gene_a and gene_b,
  # conditional on cis-SNPs for both genes
  x_a_imputed <- predict(model_a, newx = z_a, s = "lambda.1se", type = "response")
  x_b_imputed <- predict(model_b, newx = z_b, s = "lambda.1se", type = "response")
  coex <- scale((x_a - x_a_imputed) * (x_b - x_b_imputed))
  
  # Fit an elastic net model for the co-expression of gene_a and gene_b
  model_co <- cv.glmnet(x = z_both, y = coex,
                        family = "gaussian", type.measure = "mse",
                        alpha = opt$elnet_alpha, nfolds = 5, intercept = FALSE)
  
  return(list(model_a, model_b, model_co))
}


# Basic setup ------------------------------------------------------------------------------------

if (anyNA(opt[-6])) {
  stop("Some required parameters are missing. Run `cowas.R --help` for usage info.")
}

if (!dir.exists(opt$out)) {
  system(paste("mkdir", opt$out))
}

if (opt$verbose) {
  message("Processing genes ", opt$gene_a, " and ", opt$gene_b)
  message("TWAS weights and performance metrics for significant models will be saved in ", opt$out)
}


# Load and preprocess data -----------------------------------------------------------------------

if (opt$verbose) message("Importing data")

genotypes_a <- fread(file = opt$genotypes_a, sep = "\t", header = TRUE, key = "IID")
genotypes_b <- fread(file = opt$genotypes_b, sep = "\t", header = TRUE, key = "IID")
expression <- fread(file = opt$expression, sep = "\t", header = TRUE,
                    select = c("IID", opt$gene_a, opt$gene_b), key = "IID")

# Remove allele codes after each rsid, in case files were created with `plink --recode A`
setnames(genotypes_a, gsub(pattern = "_.*", replacement = "", x = names(genotypes_a)))
setnames(genotypes_b, gsub(pattern = "_.*", replacement = "", x = names(genotypes_b)))

if (is.na(opt$covariates)) {
  if (opt$verbose) message("No covariates provided. Expression will not be residualized.")
  
  # Subset to the common set of individuals
  individuals <- Reduce(intersect,
                        list(genotypes_a[, IID],
                             genotypes_b[, IID],
                             expression[, IID]))
  genotypes_a <- genotypes_a[individuals]
  genotypes_b <- genotypes_b[individuals]
  expression <- expression[individuals]
  
  # Normalize expression levels for each gene
  for (gene in c(opt$gene_a, opt$gene_b)) {
    expression_vector <- expression[[gene]]
    expression[, (gene) := scale(expression_vector)]
  }
  
} else {
  if (opt$verbose) message("Residualizing expression over covariates")
  
  covariates <- fread(file = opt$covariates, sep = "\t", header = TRUE, key = "IID")
  
  # Remove constant covariates to avoid any problems later on
  remove <- covariates[, .SD, .SDcols = !"IID"
                       ][, names(.SD), .SDcols = function(x) anyNA(x) || var(x) <= 0]
  suppressWarnings(covariates[, (remove) := NULL])
  
  # Subset to the common set of individuals (including covariates in the intersection)
  individuals <- Reduce(intersect,
                        list(genotypes_a[, IID],
                             genotypes_b[, IID],
                             expression[, IID],
                             covariates[, IID]))
  genotypes_a <- genotypes_a[individuals]
  genotypes_b <- genotypes_b[individuals]
  expression <- expression[individuals]
  covariates <- covariates[individuals]
  
  for (gene in c(opt$gene_a, opt$gene_b)) {
    regression <- covariates[, summary(lm(expression[[gene]] ~ ., data = .SD)),
                             .SDcols = !"IID"]
    expression[, (gene) := scale(regression$residuals)]
    
    if (opt$verbose) {
      message(regression$r.squared, " of variance in ", gene, " explained by covariates")
    }
  }
}

# Remove genes with no remaining variation in expression (after residualizing and normalizing)
remove <- expression[, .SD, .SDcols = !"IID"
                     ][, names(.SD), .SDcols = function(x) anyNA(x) || var(x) <= 0]
suppressWarnings(expression[, (remove) := NULL])

if (!(opt$gene_a %in% names(expression) && opt$gene_b %in% names(expression))) {
  stop("Expression of ", opt$gene_a, " or ", opt$gene_b, " is either constant or NA.
       Skipping this pair.")
}

expression_a <- expression[[opt$gene_a]]
expression_b <- expression[[opt$gene_b]]

# Subset genotypes to those present in the snp annotation matrix
snps <- fread(file = opt$snps, sep = "\t", header = FALSE)
rsids_a <- names(genotypes_a)[-1]
rsids_b <- names(genotypes_b)[-1]
snps_keep <- intersect(snps[, V2], union(rsids_a, rsids_b))
diff_a <- setdiff(rsids_a, snps_keep)
diff_b <- setdiff(rsids_b, snps_keep)
suppressWarnings(genotypes_a[, (diff_a) := NULL])
suppressWarnings(genotypes_b[, (diff_b) := NULL])

# Only keep one copy of duplicated snps, which can arise when cis-regions overlap
snps <- snps[!duplicated(V2)]

# Normalize genotype data
genotypes_a[, (rsids_a) := lapply(.SD, scale), .SDcols = rsids_a]
genotypes_b[, (rsids_b) := lapply(.SD, scale), .SDcols = rsids_b]

# Remove monomorphic snps
remove <- genotypes_a[, .SD, .SDcols = !"IID"
                      ][, names(.SD), .SDcols = function(x) anyNA(x) || var(x) <= 0]
suppressWarnings(genotypes_a[, (remove) := NULL])

remove <- genotypes_b[, .SD, .SDcols = !"IID"
                      ][, names(.SD), .SDcols = function(x) anyNA(x) || var(x) <= 0]
suppressWarnings(genotypes_b[, (remove) := NULL])

# Create a data table with all cis-SNPs for both genes
common_cols <- intersect(names(genotypes_a), names(genotypes_b))
unique_cols_b <- c("IID", setdiff(names(genotypes_b), common_cols))
genotypes_both <- genotypes_a[genotypes_b[, ..unique_cols_b],
                              on = "IID", nomatch = 0]

# Update the snp annotation matrix
snps <- snps[V2 %in% names(genotypes_both)]


# Perform CV to estimate performance -------------------------------------------------------------

# Total number of samples used
N.tot <- expression[, .N]

# Matrix for saving performance metrics
cv.performance <- matrix(data = NA, nrow = 2, ncol = 3)
rownames(cv.performance) <- c("rsq", "pval")
colnames(cv.performance) <- c("gene_a", "gene_b", "coexpression")

if (opt$cv_folds > 1) {
  
  if (opt$verbose) message("Starting cross-validation")
  
  # Shuffle the data to break correlations before defining folds
  shuffled <- sample(N.tot)
  genotypes_a_shuffled <- genotypes_a[shuffled]
  genotypes_b_shuffled <- genotypes_b[shuffled]
  genotypes_both_shuffled <- genotypes_both[shuffled]
  expression_shuffled <- expression[shuffled]
  
  folds <- cut(seq(1, N.tot), breaks = opt$cv_folds, labels = FALSE)
  
  # Matrix for storing imputed expression and co-expression in left-out folds
  imputed <- matrix(data = NA, nrow = N.tot, ncol = 3)
  colnames(imputed) <- c("gene_a", "gene_b", "coexpression")
  rownames(imputed) <- expression_shuffled[, IID]
  
  for (k in 1:opt$cv_folds) {
    test_indices <- which(folds == k, arr.ind = TRUE)
    
    genotypes_a_train <- genotypes_a_shuffled[!test_indices]
    genotypes_b_train <- genotypes_b_shuffled[!test_indices]
    genotypes_both_train <- genotypes_both_shuffled[!test_indices]
    
    genotypes_a_test <- genotypes_a_shuffled[test_indices]
    genotypes_b_test <- genotypes_b_shuffled[test_indices]
    genotypes_both_test <- genotypes_both_shuffled[test_indices]
    
    expression_a_train <- scale(expression_shuffled[!test_indices][[opt$gene_a]])
    expression_b_train <- scale(expression_shuffled[!test_indices][[opt$gene_b]])
    
    cv_cowas_models <- TrainCowasModels(genotypes_a_train,
                                        genotypes_b_train,
                                        genotypes_both_train,
                                        expression_a_train,
                                        expression_b_train)
    
    imputed[test_indices, "gene_a"] <- predict(cv_cowas_models$model_a,
                                               newx = genotypes_a_test,
                                               s = "lambda.1se",
                                               type = "response")
    imputed[test_indices, "gene_b"] <- predict(cv_cowas_models$model_b,
                                               newx = genotypes_b_test,
                                               s = "lambda.1se",
                                               type = "response")
    imputed[test_indices, "coexpression"] <- predict(cv_cowas_models$model_co,
                                                     newx = genotypes_both_test,
                                                     s = "lambda.1se",
                                                     type = "response")
    
    if (opt$verbose) message("Finished fold ", k, " of ", opt$cv_folds)
  }
  
  # Save performance metrics as long as all three models were able to predict some variation
  if (sd(imputed[, "gene_a"]) > 0 && sd(imputed[, "gene_b"]) > 0 && sd(imputed[, "coexpression"]) > 0) {
    imputed <- imputed[expression[, IID], ]
    
    a_fit <- summary(lm(expression_a ~ imputed[["gene_a"]]))
    b_fit <- summary(lm(expression_b ~ imputed[["gene_b"]]))
    
    coexpression <- scale((expression_a - imputed[["gene_a"]]) * (expression_b - imputed[["gene_b"]]))
    coexpression_fit <- summary(lm(coexpression[, 1] ~ imputed[["coexpression"]]))
    
    # R^2 of the expression and co-expression models
    cv.performance["rsq", "model_a"] <- a_fit$adj.r.squared
    cv.performance["rsq", "model_b"] <- b_fit$adj.r.squared
    cv.performance["rsq", "coexpression"] <- coexpression_fit$adj.r.squared
    
    # P-value from an F-test comparing each fitted model to a null model
    cv.performance["pval", "model_a"] <- a_fit$coefficients[2, 4]
    cv.performance["pval", "model_b"] <- b_fit$coefficients[2, 4]
    cv.performance["pval", "coexpression"] <- coexpression_fit$coefficients[2, 4]
    
    if (opt$verbose) {
      message("Gene A cross-validated R^2: ", cv.performance["rsq", "gene_a"])
      message("Gene A cross-validated F-test P-value: ", cv.performance["pval", "gene_a"])
      message("Gene B cross-validated R^2: ", cv.performance["rsq", "gene_b"])
      message("Gene B cross-validated F-test P-value: ", cv.performance["pval", "gene_b"])
      message("Co-expression cross-validated R^2: ", cv.performance["rsq", "coexpression"])
      message("Co-expression cross-validated F-test P-value: ", cv.performance["pval", "coexpression"])
    }
  } else {
    cv.performance["rsq", ] <- c(0, 0, 0)
    cv.performance["pval", ] <- c(1, 1, 1)
    
    if (opt$verbose) {
      message("All cross-validated predictions are constant. Models for this gene pair will not be saved.")
    }
  }
}


# If significant, train elastic net models on the full sample ------------------------------------

if (cv.performance["pval", "model_a"] < opt$p_threshold
    && cv.performance["pval", "model_b"] < opt$p_threshold
    && cv.performance["pval", "coexpression"] < opt$p_threshold) {
  if (opt$verbose) message("Learning full-sample eQTL weights")
  
  cowas_models <- TrainCowasModels(genotypes_a,
                                   genotypes_b,
                                   genotypes_both,
                                   expression_a,
                                   expression_b)
  
  # Don't save the intercept term because it is forced to be zero
  wgt.matrix.a <- as.matrix(coef(cowas_models$model_a, s = "lambda.1se")[-1, ])
  wgt.matrix.b <- as.matrix(coef(cowas_models$model_b, s = "lambda.1se")[-1, ])
  wgt.matrix.co <- as.matrix(coef(cowas_models$model_co, s = "lambda.1se")[-1, ])
  
  # As long as the models are not null, save all variables needed for COWAS
  if (sd(wgt.matrix.a[, 1]) > 0 && sd(wgt.matrix.b[, 1]) > 0 && sd(wgt.matrix.co[, 1]) > 0) {
    colnames(wgt.matrix.a) <- "gene_a"
    colnames(wgt.matrix.b) <- "gene_b"
    colnames(wgt.matrix.co) <- "coexpression"
    
    setDF(snps)
    save(wgt.matrix.a, wgt.matrix.b, wgt.matrix.co, snps, cv.performance, N.tot,
         file = paste0(opt$out, "/", opt$gene_a, "_", opt$gene_b, "_wgt.RData"))
  } else {
    if (opt$verbose) message("All weights are zero. This model will not be saved.")
  }
}

if (opt$verbose) message("Done!")

