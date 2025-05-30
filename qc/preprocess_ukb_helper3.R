library(data.table)

# Load UKB variants
ukb_variants <- fread(file = "data_cleaned/genotypes.pvar", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                      select = c("ID", "REF", "ALT"))

# LDL cholesterol GWAS ----------------------------------------------------------------------------

# Load GWAS data
ldl_gwas <- fread(file = "data_raw/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results",
                  header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                  select = c("rsID", "ALT", "REF", "EFFECT_SIZE", "SE", "N"))

setnames(ldl_gwas, c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n_samples"))

# Check for duplicated variants and remove any that are present in the GWAS
ldl_duplicated_variants <- ldl_gwas[duplicated(variant_id), ]$variant_id
message("\nNOTE: ", nrow(ldl_gwas[variant_id %in% ldl_duplicated_variants, ]), " rows with duplicated rsIDs removed from the Graham et al. GWAS.")
ldl_gwas <- ldl_gwas[!(variant_id %in% ldl_duplicated_variants), ]

# Subset to the variants included in both the LDL GWAS and the UKB imputed genotype data
ldl_gwas <- merge(ldl_gwas, ukb_variants, by.x = "variant_id", by.y = "ID", sort = FALSE)
ldl_gwas <- ldl_gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]

# Compute Z-scores because they were not provided
ldl_gwas[, z_score := as.numeric(beta) / as.numeric(standard_error)]
ldl_gwas[, c("beta", "standard_error") := NULL]

# Flip the LDL GWAS summary statistics when necessary
ldl_gwas[effect_allele != ALT, flip := TRUE]
ldl_gwas[flip == TRUE, `:=` (effect_allele_flipped = other_allele,
                             other_allele_flipped = effect_allele,
                             z_score = -1 * as.numeric(z_score))]
ldl_gwas[flip == TRUE, `:=` (effect_allele = effect_allele_flipped,
                             other_allele = other_allele_flipped)]
ldl_gwas[, c("REF", "ALT", "flip", "effect_allele_flipped", "other_allele_flipped") := NULL]

# The final list of variants that will be used for LDL analysis
ldl_gwas <- na.omit(ldl_gwas)
ldl_variants_keep <- ldl_gwas$variant_id

# Save the processed GWAS summary data
write.table(ldl_variants_keep, file = "TEMP_mutual_LDL_variants.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(ldl_gwas, file = "data_cleaned/Graham_2021_LDL_GWAS.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

rm(ldl_gwas, ldl_duplicated_variants, ldl_variants_keep)

# EADB GWAS of Alzheimer's disease ----------------------------------------------------------------

# Load GWAS data
eadb_gwas <- fread(file = "data_raw/GCST90027158_buildGRCh38.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                 select = c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n_cases", "n_controls"))

# Check for duplicated variants and remove any that are present in the GWAS
eadb_duplicated_variants <- eadb_gwas[duplicated(variant_id), ]$variant_id
message("NOTE: ", nrow(eadb_gwas[variant_id %in% eadb_duplicated_variants, ]), " rows with duplicated rsIDs removed from the Bellenguez et al. GWAS.")
eadb_gwas <- eadb_gwas[!(variant_id %in% eadb_duplicated_variants), ]

# Subset to the variants included in both the EADB GWAS and the UKB imputed genotype data
eadb_gwas <- merge(eadb_gwas, ukb_variants, by.x = "variant_id", by.y = "ID", sort = FALSE)
eadb_gwas <- eadb_gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]

# Compute Z-scores because they were not provided
eadb_gwas[, z_score := as.numeric(beta) / as.numeric(standard_error)]
eadb_gwas[, c("beta", "standard_error") := NULL]

# Flip the EADB GWAS summary statistics when necessary
eadb_gwas[effect_allele != ALT, flip := TRUE]
eadb_gwas[flip == TRUE, `:=` (effect_allele_flipped = other_allele,
                              other_allele_flipped = effect_allele,
                              z_score = -1 * as.numeric(z_score))]
eadb_gwas[flip == TRUE, `:=` (effect_allele = effect_allele_flipped,
                              other_allele = other_allele_flipped)]
eadb_gwas[, c("REF", "ALT", "flip", "effect_allele_flipped", "other_allele_flipped") := NULL]

# Create a column for the total sample size
eadb_gwas[, n_samples := as.numeric(n_cases) + as.numeric(n_controls)]
eadb_gwas[, c("n_cases", "n_controls") := NULL]

# This is the final list of variants that will be used for COWAS
eadb_gwas <- na.omit(eadb_gwas)
eadb_variants_keep <- eadb_gwas$variant_id

# Save the processed EADB GWAS summary data
write.table(eadb_variants_keep, file = "TEMP_mutual_AD_EADB_variants.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(eadb_gwas, file = "data_cleaned/Bellenguez_2022_AD_GWAS.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

rm(eadb_gwas, eadb_duplicated_variants, eadb_variants_keep)

# IGAP GWAS of Alzheimer's disease ----------------------------------------------------------------

# Load IGAP Alzheimer's disease GWAS data
igap_gwas <- fread(file = "data_raw/Kunkle_etal_Stage1_results.txt", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                   select = c("MarkerName", "Effect_allele", "Non_Effect_allele", "Beta", "SE"))
setnames(igap_gwas, c("variant_id", "effect_allele", "other_allele", "beta", "standard_error"))

# Check for duplicated variants and remove any that are present in the GWAS
igap_duplicated_variants <- igap_gwas[duplicated(variant_id), ]$variant_id
message("NOTE: ", nrow(igap_gwas[variant_id %in% igap_duplicated_variants, ]), " rows with duplicated rsIDs removed from the Kunkle et al. GWAS.")
igap_gwas <- igap_gwas[!(variant_id %in% igap_duplicated_variants), ]

# Subset to the variants included in both the IGAP GWAS and the UKB imputed genotype data
igap_gwas <- merge(igap_gwas, ukb_variants, by.x = "variant_id", by.y = "ID", sort = FALSE)
igap_gwas <- igap_gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]

# Compute Z-scores because they were not provided
igap_gwas[, z_score := as.numeric(beta) / as.numeric(standard_error)]
igap_gwas[, c("beta", "standard_error") := NULL]

# Flip the IGAP GWAS summary statistics when necessary
igap_gwas[effect_allele != ALT, flip := TRUE]
igap_gwas[flip == TRUE, `:=` (effect_allele_flipped = other_allele,
                              other_allele_flipped = effect_allele,
                              z_score = -1 * as.numeric(z_score))]
igap_gwas[flip == TRUE, `:=` (effect_allele = effect_allele_flipped,
                              other_allele = other_allele_flipped)]
igap_gwas[, c("REF", "ALT", "flip", "effect_allele_flipped", "other_allele_flipped") := NULL]

# Create a column for the total sample size
igap_gwas[, n_samples := 63926]

# This is the final list of variants that will be used for COWAS
igap_gwas <- na.omit(igap_gwas)
igap_variants_keep <- igap_gwas$variant_id

# Save the processed IGAP GWAS summary data
write.table(igap_variants_keep, file = "TEMP_mutual_AD_IGAP_variants.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(igap_gwas, file = "data_cleaned/Kunkle_2019_AD_GWAS.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

rm(igap_gwas, igap_duplicated_variants, igap_variants_keep)

# Parkinson's disease GWAS ------------------------------------------------------------------------

# Load GWAS data
pd_gwas <- fread(file = "data_raw/GCST009325.h.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                 select = c("rsid", "effect_allele", "other_allele", "beta", "standard_error", "N_cases", "N_controls"))

setnames(pd_gwas, c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n_cases", "n_controls"))

# Check for duplicated variants and remove any that are present in the GWAS
pd_duplicated_variants <- pd_gwas[duplicated(variant_id), ]$variant_id
message("NOTE: ", nrow(pd_gwas[variant_id %in% pd_duplicated_variants, ]), " rows with duplicated rsIDs removed from the Nalls et al. GWAS.\n")
pd_gwas <- pd_gwas[!(variant_id %in% pd_duplicated_variants), ]

# Subset to the variants included in both the PD GWAS and the UKB imputed genotype data
pd_gwas <- merge(pd_gwas, ukb_variants, by.x = "variant_id", by.y = "ID", sort = FALSE)
pd_gwas <- pd_gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]

# Compute Z-scores because they were not provided
pd_gwas[, z_score := as.numeric(beta) / as.numeric(standard_error)]
pd_gwas[, c("beta", "standard_error") := NULL]

# Flip the PD GWAS summary statistics when necessary
pd_gwas[effect_allele != ALT, flip := TRUE]
pd_gwas[flip == TRUE, `:=` (effect_allele_flipped = other_allele,
                            other_allele_flipped = effect_allele,
                            z_score = -1 * as.numeric(z_score))]
pd_gwas[flip == TRUE, `:=` (effect_allele = effect_allele_flipped,
                            other_allele = other_allele_flipped)]
pd_gwas[, c("REF", "ALT", "flip", "effect_allele_flipped", "other_allele_flipped") := NULL]

# Create a column for the total sample size
pd_gwas[, n_samples := as.numeric(n_cases) + as.numeric(n_controls)]
pd_gwas[, c("n_cases", "n_controls") := NULL]

# The final list of variants that will be used for PD analysis
pd_gwas <- na.omit(pd_gwas)
pd_variants_keep <- pd_gwas$variant_id

# Save the processed GWAS summary data
write.table(pd_variants_keep, file = "TEMP_mutual_PD_variants.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(pd_gwas, file = "data_cleaned/Nalls_2019_PD_GWAS.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

rm(pd_gwas, pd_duplicated_variants, pd_variants_keep)

