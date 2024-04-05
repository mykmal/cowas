library(data.table)

# Load UKB variants
ukb_variants <- fread(file = "data_cleaned/genotypes.pvar", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                      select = c("ID", "REF", "ALT"))

# Alzheimer's disease GWAS ------------------------------------------------------------------------

# Load GWAS data
ad_gwas <- fread(file = "data_raw/GCST90027158_buildGRCh38.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                 select = c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n_cases", "n_controls"))

# The Bellenguez et al. GWAS has some rows with identical rsIDs, positions, effect alleles, and alternate alleles but different p-values and allele frequencies.
# This seems to be a mistake in the data, so we remove all such rows.
ad_duplicated_variants <- ad_gwas[duplicated(variant_id), ]$variant_id
message("\nNOTE: ", nrow(ad_gwas[variant_id %in% ad_duplicated_variants, ]), " duplicated rows removed from the Bellenguez et al. GWAS\n")
ad_gwas <- ad_gwas[!variant_id %in% ad_duplicated_variants, ]

# Subset to the variants included in both the AD GWAS and the UKB imputed genotype data
ad_gwas <- merge(ad_gwas, ukb_variants, by.x = "variant_id", by.y = "ID", sort = FALSE)
ad_gwas <- ad_gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]

# Compute Z-scores because they were not provided by Bellenguez et al.
ad_gwas[, z_score := as.numeric(beta) / as.numeric(standard_error)]
ad_gwas[, c("beta", "standard_error") := NULL]

# Flip the AD GWAS summary statistics when necessary
ad_gwas[effect_allele != ALT, flip := TRUE]
ad_gwas[flip == TRUE, `:=` (effect_allele_flipped = other_allele,
                            other_allele_flipped = effect_allele,
                            z_score = -1 * as.numeric(z_score))]
ad_gwas[flip == TRUE, `:=` (effect_allele = effect_allele_flipped,
                            other_allele = other_allele_flipped)]
ad_gwas[, c("REF", "ALT", "flip", "effect_allele_flipped", "other_allele_flipped") := NULL]

# Create a column for the total sample size
ad_gwas[, n_samples := as.numeric(n_cases) + as.numeric(n_controls)]
ad_gwas[, c("n_cases", "n_controls") := NULL]

# The final list of variants that will be used for AD analysis
ad_gwas <- na.omit(ad_gwas)
ad_variants_keep <- ad_gwas$variant_id

# Save the processed GWAS summary data
write.table(ad_variants_keep, file = "TEMP_mutual_AD_variants.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(ad_gwas, file = "data_cleaned/Bellenguez_2022_AD_GWAS.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

rm(ad_gwas, ad_duplicated_variants, ad_variants_keep)

# Schizophrenia GWAS ------------------------------------------------------------------------------

# Load GWAS data
scz_gwas <- fread(file = "data_raw/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv",
                  skip = 73, header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                  select = c("ID", "A1", "A2", "BETA", "SE", "NCAS", "NCON"))

setnames(scz_gwas, c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n_cases", "n_controls"))

# Subset to the variants included in both the SCZ GWAS and the UKB imputed genotype data
scz_gwas <- merge(scz_gwas, ukb_variants, by.x = "variant_id", by.y = "ID", sort = FALSE)
scz_gwas <- scz_gwas[(effect_allele == ALT & other_allele == REF) | (effect_allele == REF & other_allele == ALT), ]

# Compute Z-scores because they were not provided by Trubetskoy et al.
scz_gwas[, z_score := as.numeric(beta) / as.numeric(standard_error)]
scz_gwas[, c("beta", "standard_error") := NULL]

# Flip the SCZ GWAS summary statistics when necessary
scz_gwas[effect_allele != ALT, flip := TRUE]
scz_gwas[flip == TRUE, `:=` (effect_allele_flipped = other_allele,
                             other_allele_flipped = effect_allele,
                             z_score = -1 * as.numeric(z_score))]
scz_gwas[flip == TRUE, `:=` (effect_allele = effect_allele_flipped,
                             other_allele = other_allele_flipped)]
scz_gwas[, c("REF", "ALT", "flip", "effect_allele_flipped", "other_allele_flipped") := NULL]

# Create a column for the total sample size
scz_gwas[, n_samples := as.numeric(n_cases) + as.numeric(n_controls)]
scz_gwas[, c("n_cases", "n_controls") := NULL]

# The final list of variants that will be used for SCZ analysis
scz_gwas <- na.omit(scz_gwas)
scz_variants_keep <- scz_gwas$variant_id

# Save the processed GWAS summary data
write.table(scz_variants_keep, file = "TEMP_mutual_SCZ_variants.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(scz_gwas, file = "data_cleaned/Trubetskoy_2022_SCZ_GWAS.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

