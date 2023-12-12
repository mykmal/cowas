library(data.table)

ukb_variants <- fread(file = "TEMP_all_ukb_variants.pvar", header = TRUE, stringsAsFactors = FALSE, select = c("ID", "REF", "ALT"))
gwas <- fread(file = "raw/GCST90027158_buildGRCh38.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)

gwas_variants <- gwas[, c("variant_id", "effect_allele", "other_allele")]

mutual_variants <- merge(ukb_variants, gwas_variants, by.x = "ID", by.y = "variant_id", sort = FALSE)
mutual_variants <- mutual_variants[(REF == effect_allele & ALT == other_allele) | (REF == other_allele & ALT == effect_allele)]
mutual_variants <- mutual_variants$ID

gwas <- gwas[variant_id %in% mutual_variants]

fwrite(mutual_variants, file = "TEMP_mutual_variants.txt", quote = FALSE, col.names = FALSE, eol = "\n", compress = "none")
fwrite(gwas, file = "gwas/Bellenguez_2022_AD_gwas.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

