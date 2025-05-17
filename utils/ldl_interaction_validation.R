# This script uses individual-level data from the UK Biobank to regress LDL cholesterol levels on protein expression interactions.
# In particular, we fit the model LDL ~ A + B + A * B for each well-imputed pair A, B.
# We performed this analysis to validate COWAS results for LDL cholesterol.

library(data.table)

# Load the cleaned COWAS association results for LDL, as processed by create_plots.Rmd
cowas_results <- fread(file = "figures_tables/csv/results_ldl.csv", header = TRUE, sep = ",", na.strings = "NA", stringsAsFactors = FALSE)
n_pairs <- nrow(cowas_results)

# Create an empty data frame for storing association results
results <- data.frame("PROTEIN_A" = character(n_pairs),
                      "PROTEIN_B" = character(n_pairs),
                      "SAMPLE_SIZE" = numeric(n_pairs),
                      "BETA_1" = numeric(n_pairs),
                      "SE_BETA_1" = numeric(n_pairs),
                      "PVAL_BETA_1" = numeric(n_pairs),
                      "BETA_2" = numeric(n_pairs),
                      "SE_BETA_2" = numeric(n_pairs),
                      "PVAL_BETA_2" = numeric(n_pairs),
                      "BETA_3" = numeric(n_pairs),
                      "SE_BETA_3" = numeric(n_pairs),
                      "PVAL_BETA_3" = numeric(n_pairs))

# Load LDL cholesterol measurements from the UK Biobank phase 2 metabolomics dataset (data field 23405)
ldl <- fread(file = "data_raw/metabolites_phase2.tab", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
             select = c("f.eid", "f.23405.0.0"))
ldl <- na.omit(ldl)
setnames(ldl, c("IID", "LDL"))

# Load protein expression data
expression <- fread(file = "data_cleaned/proteins.tsv", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)

# Load covariate data
covariates <- fread(file = "data_cleaned/covariates.tsv", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
covariates <- na.omit(covariates)

# Loop through the list of well-imputed protein pairs
for (pair_number in 1:n_pairs) {
  
  protein_a <- cowas_results[pair_number, ]$ID_A
  protein_b <- cowas_results[pair_number, ]$ID_B
  
  # Subset the expression data table to the current pair of proteins
  current_expression <- subset(expression, select = c("IID", protein_a, protein_b))
  current_expression <- na.omit(current_expression)
  
  # Subset all three datasets to the common set of individuals with no missing values for the current pair
  current_data <- merge(ldl, current_expression, by = "IID")
  current_data <- merge(current_data, covariates, by = "IID")
  
  # Temporarily replace hyphens in protein names with underscores, because R does not like hyphens
  setnames(current_data, gsub("-", "_", names(current_data), fixed = TRUE))
  
  # Center and scale each covariate
  for (column in names(covariates)[-1]) {
    set(x = current_data, j = column, value = scale(current_data[[column]]))
  }
  
  # Process each of the two proteins
  for (protein in c(protein_a, protein_b)) {
    
    # Temporarily replace hyphens in protein names with underscores, because R does not like hyphens
    protein <- gsub("-", "_", protein, fixed = TRUE)
    
    # Compute the rank of each observation
    ranks <- rank(current_data[[protein]], ties.method = "average")
    
    # Perform a rank-based inverse normal transformation
    current_data[, (protein) := stats::qnorm((ranks - 0.375) / (nrow(current_data) - 2 * 0.375 + 1))]
    
    # Regress out covariates
    formula <- as.formula(paste0(protein, " ~ ", paste0(names(covariates)[-1], collapse = " + ")))
    regression <- stats::lm(formula, data = current_data)
    current_data[, (protein) := regression$residuals]
  }
  
  # Perform a rank-based inverse normal transformation for LDL cholesterol
  ranks <- rank(current_data$LDL, ties.method = "average")
  current_data[, LDL := stats::qnorm((ranks - 0.375) / (nrow(current_data) - 2 * 0.375 + 1))]
  
  # Regress out covariates for LDL cholesterol
  formula <- as.formula(paste0("LDL ~ ", paste0(names(covariates)[-1], collapse = " + ")))
  regression <- stats::lm(formula, data = current_data)
  current_data[, LDL := regression$residuals]
  
  # Finally, regress LDL residuals on the interaction of protein expression residuals
  # Namely, we fit the model LDL ~ A + B + A*B
  formula <- as.formula(paste0("LDL ~ ", gsub("-", "_", protein_a, fixed = TRUE),
                               " * ", gsub("-", "_", protein_b, fixed = TRUE)))
  model <- summary(stats::lm(formula, data = current_data))
  
  # Save the association statistics
  results[pair_number, 1:2] <- c(protein_a, protein_b)
  results[pair_number, 3:12] <- c(nrow(current_data),
                                  model$coefficients[2, 1],
                                  model$coefficients[2, 2],
                                  model$coefficients[2, 4],
                                  model$coefficients[3, 1],
                                  model$coefficients[3, 2],
                                  model$coefficients[3, 4],
                                  model$coefficients[4, 1],
                                  model$coefficients[4, 2],
                                  model$coefficients[4, 4])
}

# Compute the percentage of significant COWAS interaction tests that were validated by this analysis
sig_cowas_interaction <- cowas_results[PVAL_THETA_JOINT_CO < 0.05 / n_pairs, ]
sig_cowas_interaction[, ID := paste0(ID_A, "_", ID_B)]

results$ID <- paste0(results$PROTEIN_A, "_", results$PROTEIN_B)
sig_validation <- subset(results, PVAL_BETA_3 < 0.05 & ID %in% sig_cowas_interaction$ID)
print(nrow(sig_validation) / nrow(sig_cowas_interaction))

# Save the results
results <- results[, -13]
write.csv(results, file = "association_results/ldl_validation.csv", quote = FALSE, row.names = FALSE)

