library(data.table)

# Load data
ukb_variants <- fread(file = "TEMP_merged.pvar", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                      select = c("ID", "REF", "ALT"))
gwas <- fread(file = "data_raw/GCST90027158_buildGRCh38.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
              select = c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n_cases", "n_controls"))

# The Bellenguez et al. GWAS has some rows with identical rsIDs, positions, effect alleles, and alternate alleles but different p-values and allele frequencies
# This seems to be a mistake in the data, so we remove all such rows
duplicated_variants <- gwas[duplicated(gwas[, c("variant_id", "effect_allele", "other_allele")]), ]$variant_id
message("\nNOTE: ", nrow(gwas[variant_id %in% duplicated_variants, ]), " duplicated rows removed from the Bellenguez et al. GWAS\n")
gwas <- gwas[!variant_id %in% duplicated_variants, ]

# Subset to the variants included in both the GWAS and the UKB imputed genotype data
gwas_mutual_variants <- merge(gwas, ukb_variants, by.x = "variant_id", by.y = "ID", sort = FALSE)
gwas_mutual_variants <- gwas_mutual_variants[(effect_allele == REF & other_allele == ALT) | (effect_allele == ALT & other_allele == REF), ]
gwas_mutual_variants[, c("REF", "ALT") := NULL]

# Compute Z-scores because they were not provided by Bellenguez et al.
gwas_mutual_variants[, z_score := as.numeric(beta) / as.numeric(standard_error)]
gwas_mutual_variants[, c("beta", "standard_error") := NULL]

# The final list of variants that will be used for analysis
variants_keep <- gwas_mutual_variants$variant_id

# Save the processed GWAS summary data
write.table(variants_keep, file = "TEMP_mutual_variants.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(gwas_mutual_variants, file = "data_cleaned/Bellenguez_2022_AD_gwas.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

