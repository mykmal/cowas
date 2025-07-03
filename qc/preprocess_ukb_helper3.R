library(data.table)

# Define the steps for processing GWAS summary statistics
process_gwas <- function(filename, column_names, label) {
  
  # Load the GWAS data
  gwas <- fread(file = paste0("data_raw/", filename),
                header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                select = column_names)
  
  if (ncol(gwas) == 5) {
    # In this case, we are dealing with the AD IGAP GWAS
    
    # Standardize the column names
    setnames(gwas, c("variant_id", "effect_allele", "other_allele", "beta", "standard_error"))
    
    # Create a column for the total sample size
    gwas[, n_samples := 63926]
    
  } else if (ncol(gwas) == 6) {
    # In this case, we are dealing with the LDL GWAS
    
    # Standardize the column names
    setnames(gwas, c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n_samples"))
    
  } else if (ncol(gwas) == 7) {
    # In this case, we are dealing with the AD EADB GWAS or the PD GWAS
    
    # Standardize the column names
    setnames(gwas, c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n_cases", "n_controls"))
    
    # Create a column for the total sample size
    gwas[, n_samples := as.numeric(n_cases) + as.numeric(n_controls)]
    gwas[, c("n_cases", "n_controls") := NULL] 
  }
  
  # Check for duplicated variants and remove any that are present in the GWAS
  duplicated_variants <- gwas[duplicated(variant_id), ]$variant_id
  message("NOTE: ", nrow(gwas[variant_id %in% duplicated_variants, ]), " rows with duplicated rsIDs removed from the ", label, " GWAS.")
  gwas <- gwas[!(variant_id %in% duplicated_variants), ]
  
  # Subset to the variants included in both the GWAS and in the UKB imputed genotype data
  gwas <- merge(gwas, ukb_variants, by.x = "variant_id", by.y = "ID", sort = FALSE)
  gwas <- gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]
  
  # Compute Z scores because they were not provided
  gwas[, z_score := as.numeric(beta) / as.numeric(standard_error)]
  gwas[, c("beta", "standard_error") := NULL]
  
  # Flip the GWAS summary statistics when necessary
  gwas[effect_allele != ALT, flip := TRUE]
  gwas[flip == TRUE, `:=` (effect_allele_flipped = other_allele,
                           other_allele_flipped = effect_allele,
                           z_score = -1 * as.numeric(z_score))]
  gwas[flip == TRUE, `:=` (effect_allele = effect_allele_flipped,
                           other_allele = other_allele_flipped)]
  gwas[, c("REF", "ALT", "flip", "effect_allele_flipped", "other_allele_flipped") := NULL]
  
  # This is the final list of variants that will be used
  gwas <- na.omit(gwas)
  variants_keep <- gwas$variant_id
  
  # Save the processed GWAS summary data
  write.table(variants_keep, file = paste0("TEMP_mutual_", label, "_variants.txt"),
              quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
  fwrite(gwas, file = paste0("data_cleaned/", label, "_GWAS.tsv"),
         quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")
}

# Load UKB variants
ukb_variants <- fread(file = "data_cleaned/genotypes.pvar", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                      select = c("ID", "REF", "ALT"))

# Print a newline so the log file looks nice
message("\n")

# Process the summary statistics for each GWAS
process_gwas(filename = "LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results",
             column_names = c("rsID", "ALT", "REF", "EFFECT_SIZE", "SE", "N"),
             label = "LDL")

process_gwas(filename = "GCST90027158_buildGRCh38.tsv",
             column_names = c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n_cases", "n_controls"),
             label = "AD_EADB")

process_gwas(filename = "Kunkle_etal_Stage1_results.txt",
             column_names = c("MarkerName", "Effect_allele", "Non_Effect_allele", "Beta", "SE"),
             label = "AD_IGAP")

process_gwas(filename = "GCST009325.h.tsv",
             column_names = c("rsid", "effect_allele", "other_allele", "beta", "standard_error", "N_cases", "N_controls"),
             label = "PD")

# Print a newline so the log file looks nice
message("\n")

