# This script calculates the average physical distance among proteins identified as significant
# by the COWAS global or interaction tests, and the average physical distance among those that
# were not significant. It also calculates the average observed co-expression between protein
# pairs in each of these sets.

library(data.table)

# Calculate the average distance between protein pairs in the provided data table
calculate_distance <- function(results) {
  
  # Number of pairs not on the same chromosome
  n_diff_chr <- 0
  # Number of pairs on the same chromosome
  n_same_chr <- 0
  # Total distance among pairs on the same chromosome
  total_dist <- 0
  
  # Loop through all protein pairs
  for (pair in 1:nrow(results)) {
    
    protein_a <- results[pair, ]$ID_A
    protein_b <- results[pair, ]$ID_B
    
    # Get the chromosome of the gene coding each protein
    chr_a <- annotations[Assay == protein_a, ][1, ]$chr_hg19
    chr_b <- annotations[Assay == protein_b, ][1, ]$chr_hg19
    
    # If the proteins are coded by genes on different chromosomes, we can't calculate a distance.
    # Record this and move on to the next pair.
    if (chr_a != chr_b) {
      n_diff_chr <- n_diff_chr + 1
      next
    }
    
    # Get the start position of the gene coding each protein
    start_a <- annotations[Assay == protein_a, ][1, ]$start_hg19
    start_b <- annotations[Assay == protein_b, ][1, ]$start_hg19
    
    # Calculate the physical distance between these genes
    total_dist <- total_dist + abs(start_a - start_b)
    n_same_chr <- n_same_chr + 1
  }
  
  avg_dist <- total_dist / n_same_chr
  
  return(list(n_diff_chr = n_diff_chr,
              n_same_chr = n_same_chr,
              avg_dist = avg_dist))
}

# Calculate the average Pearson correlation between protein pairs in the provided data table
calculate_correlation <- function(results) {
  
  correlations <- 0
  
  # Loop through all protein pairs
  for (pair in 1:nrow(results)) {
    
    protein_a <- results[pair, ]$ID_A
    protein_b <- results[pair, ]$ID_B
    
    # Subset expression data to samples with no missing values for the current pair
    expression_current <- subset(full_expression, select = c("IID", protein_a, protein_b))
    expression_current <- na.omit(expression_current)
    
    # Subset expression data and covariates to their common set of samples
    current_data <- merge(expression_current, covariates, by = "IID")
    
    # Replace hyphens in protein names with underscores, because R does not like hyphens
    setnames(current_data, gsub("-", "_", names(current_data), fixed = TRUE))
    
    # Center and scale each covariate
    for (column in names(covariates)[-1]) {
      set(x = current_data, j = column, value = scale(current_data[[column]]))
    }
    
    # Replace hyphens in protein names with underscores, because R does not like hyphens
    protein_a <- gsub("-", "_", protein_a, fixed = TRUE)
    protein_b <- gsub("-", "_", protein_b, fixed = TRUE)
    
    # Process each of the two proteins
    for (protein in c(protein_a, protein_b)) {
      
      # Compute the rank of each observation
      ranks <- rank(current_data[[protein]], ties.method = "average")
      
      # Perform a rank-based inverse normal transformation
      current_data[, (protein) := stats::qnorm((ranks - 0.375) / (nrow(current_data) - 2 * 0.375 + 1))]
      
      # Regress out covariates
      formula <- as.formula(paste0(protein, " ~ ", paste0(names(covariates)[-1], collapse = " + ")))
      regression <- stats::lm(formula, data = current_data)
      current_data[, (protein) := regression$residuals]
    }
    
    # Save the correlation between the two proteins in this pair
    correlations <- cbind(correlations, cor(current_data[[protein_a]],
                                            current_data[[protein_b]],
                                            method = "pearson"))
  }
  
  # Remove the zero value that we started with
  correlations <- correlations[-1]
  
  # Return the average absolute value correlation
  avg_correlation <- mean(abs(correlations))
  
  return(avg_correlation)
}

# Load protein annotations
annotations <- fread(file = "3_olink_annotations_lifted.tsv", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)

# Load cleaned COWAS association results for each trait, as processed by create_plots.Rmd
ldl_results <- fread(file = "figures_tables/csv/results_ldl.csv", header = TRUE, sep = ",", na.strings = "NA", stringsAsFactors = FALSE)
ad_results <- fread(file = "figures_tables/csv/results_ad.csv", header = TRUE, sep = ",", na.strings = "NA", stringsAsFactors = FALSE)
ad_igap_results <- fread(file = "figures_tables/csv/results_ad_igap.csv", header = TRUE, sep = ",", na.strings = "NA", stringsAsFactors = FALSE)
ad_product_results <- fread(file = "figures_tables/csv/results_ad_product.csv", header = TRUE, sep = ",", na.strings = "NA", stringsAsFactors = FALSE)
pd_results <- fread(file = "figures_tables/csv/results_pd.csv", header = TRUE, sep = ",", na.strings = "NA", stringsAsFactors = FALSE)

# Load protein expression data
full_expression <- fread(file = "data_cleaned/proteins.tsv", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)

# Load covariate data
covariates <- fread(file = "data_cleaned/covariates.tsv", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
covariates <- na.omit(covariates)

# Subset each results dataset to significant and non-significant pairs
ldl_sig <- ldl_results[PVAL_THETA_JOINT_CO < 0.05 / nrow(ldl_results) | PVAL_FSTAT_JOINT < 0.05 / nrow(ldl_results), ]
ldl_nonsig <- ldl_results[!paste0(ID_A, ID_B) %in% paste0(ldl_sig$ID_A, ldl_sig$ID_B), ]

ad_sig <- ad_results[PVAL_THETA_JOINT_CO < 0.05 / nrow(ad_results) | PVAL_FSTAT_JOINT < 0.05 / nrow(ad_results), ]
ad_nonsig <- ad_results[!paste0(ID_A, ID_B) %in% paste0(ad_sig$ID_A, ad_sig$ID_B), ]

ad_igap_sig <- ad_igap_results[PVAL_THETA_JOINT_CO < 0.05 / nrow(ad_igap_results) | PVAL_FSTAT_JOINT < 0.05 / nrow(ad_igap_results), ]
ad_igap_nonsig <- ad_igap_results[!paste0(ID_A, ID_B) %in% paste0(ad_igap_sig$ID_A, ad_igap_sig$ID_B), ]

ad_prod_sig <- ad_product_results[PVAL_THETA_JOINT_CO < 0.05 / nrow(ad_product_results) | PVAL_FSTAT_JOINT < 0.05 / nrow(ad_product_results), ]
ad_prod_nonsig <- ad_product_results[!paste0(ID_A, ID_B) %in% paste0(ad_prod_sig$ID_A, ad_prod_sig$ID_B), ]

# Remove GLIPR1 because it only has 106 samples and covariate adjustment fails for it
ad_prod_nonsig <- ad_prod_nonsig[ID_B != "GLIPR1", ]

pd_sig <- pd_results[PVAL_THETA_JOINT_CO < 0.05 / nrow(pd_results) | PVAL_FSTAT_JOINT < 0.05 / nrow(pd_results), ]
pd_nonsig <- pd_results[!paste0(ID_A, ID_B) %in% paste0(pd_sig$ID_A, pd_sig$ID_B), ]

# Calculate the average distance in each set
ldl_dist_sig <- calculate_distance(ldl_sig)
ldl_dist_nonsig <- calculate_distance(ldl_nonsig)

ad_dist_sig <- calculate_distance(ad_sig)
ad_dist_nonsig <- calculate_distance(ad_nonsig)

ad_igap_dist_sig <- calculate_distance(ad_igap_sig)
ad_igap_dist_nonsig <- calculate_distance(ad_igap_nonsig)

ad_prod_dist_sig <- calculate_distance(ad_prod_sig)
ad_prod_dist_nonsig <- calculate_distance(ad_prod_nonsig)

pd_dist_sig <- calculate_distance(pd_sig)
pd_dist_nonsig <- calculate_distance(pd_nonsig)

# Calculate the average correlation in each set
ldl_cor_sig <- calculate_correlation(ldl_sig)
ldl_cor_nonsig <- calculate_correlation(ldl_nonsig)

ad_cor_sig <- calculate_correlation(ad_sig)
ad_cor_nonsig <- calculate_correlation(ad_nonsig)

ad_igap_cor_sig <- calculate_correlation(ad_igap_sig)
ad_igap_cor_nonsig <- calculate_correlation(ad_igap_nonsig)

ad_prod_cor_sig <- calculate_correlation(ad_prod_sig)
ad_prod_cor_nonsig <- calculate_correlation(ad_prod_nonsig)

pd_cor_sig <- calculate_correlation(pd_sig)
pd_cor_nonsig <- calculate_correlation(pd_nonsig)

