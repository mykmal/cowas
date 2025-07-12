# This script extracts a set of initial predictors to be used when training COWAS imputation models.

library(data.table)

# Get the list of protein names
sumstat_files <- list.files("pqtl_associations/", pattern = "*.sumstats.tsv")
protein_names <- gsub(".sumstats.tsv", "", sumstat_files, fixed = TRUE)

dir.create("general_predictors")

# Part 1: extract predictors from among all variants

dir.create("general_predictors/predictors_top_global_beta")
dir.create("general_predictors/predictors_top_global_pval")

for (protein in protein_names) {
  sumstats <- fread(file = paste0("pqtl_associations/", protein, ".sumstats.tsv"),
                    header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
                    select = c("ID", "BETA", "P"), verbose = FALSE)
  
  # Select 100 variants with the largest absolute value betas
  sumstats[, BETA := abs(BETA)]
  sumstats <- sumstats[order(sumstats$BETA, decreasing = TRUE), ]
  sumstats_top <- sumstats[1:100, ]
  write.table(t(sumstats_top$ID), file = paste0("general_predictors/predictors_top_global_beta/", protein, ".variants.txt"),
              quote = FALSE, sep = "\n", eol = "\n", row.names = FALSE, col.names = FALSE)
  
  # Select 100 variants with the smallest P values
  sumstats[, P := as.numeric(P)]
  sumstats <- sumstats[order(sumstats$P, decreasing = FALSE), ]
  sumstats_top <- sumstats[1:100, ]
  write.table(t(sumstats_top$ID), file = paste0("general_predictors/predictors_top_global_pval/", protein, ".variants.txt"),
              quote = FALSE, sep = "\n", eol = "\n", row.names = FALSE, col.names = FALSE)
}


# Part 2: extract predictors from variants in the cis region

dir.create("general_predictors/predictors_top_cis_beta")
dir.create("general_predictors/predictors_top_cis_pval")

annotations <- fread(file = "protein_annotations_derived/3_olink_annotations_lifted.tsv",
                     header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
                     select = c("Assay", "chr_hg19", "start_hg19", "end_hg19"), verbose = FALSE)

# Exclude proteins coded by non-autosomal genes
annotations <- annotations[chr_hg19 != "X", ]
protein_names <- protein_names[protein_names %in% annotations$Assay]

# Define the cis regions of each protein-coding gene
annotations[, start_hg19 := start_hg19 - 500000]
annotations[, end_hg19 := end_hg19 + 500000]

for (protein in protein_names) {
  # Get annotations for the current protein
  annotations_subset <- annotations[Assay == protein, ]
  
  # If more than one gene matches, use the first one listed
  if (nrow(annotations_subset) > 1) {
    annotations_subset <- annotations_subset[1, ]
  }
  
  # Read pQTL summary statistics for the current protein
  sumstats <- fread(file = paste0("pqtl_associations/", protein, ".sumstats.tsv"),
                    header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
                    select = c("#CHROM", "POS", "ID", "BETA", "P"), verbose = FALSE)
  setnames(sumstats, "#CHROM", "CHROM")
  
  # Subset the pQTL summary statistics to variants in the cis region of the gene coding this protein
  sumstats <- sumstats[(CHROM == annotations_subset$chr_hg19) & (POS >= annotations_subset$start) & (POS <= annotations_subset$end), ]
  
  # Report proteins that have less than 100 cis-variants
  if (nrow(sumstats) == 0) {
    message(protein, " has no cis-variants")
    next
  } else if (nrow(sumstats) < 100) {
    n_snps <- nrow(sumstats)
    message(protein, " has only ", n_snps, " cis-variants")
  } else {
    n_snps <- 100
  }
  
  # Select up to 100 cis-variants with the largest absolute value betas
  sumstats[, BETA := abs(BETA)]
  sumstats <- sumstats[order(sumstats$BETA, decreasing = TRUE), ]
  sumstats_top <- sumstats[1:n_snps, ]
  write.table(t(sumstats_top$ID), file = paste0("general_predictors/predictors_top_cis_beta/", protein, ".variants.txt"),
              quote = FALSE, sep = "\n", eol = "\n", row.names = FALSE, col.names = FALSE)
  
  # Select up to 100 cis-variants with the smallest P values
  sumstats[, P := as.numeric(P)]
  sumstats <- sumstats[order(sumstats$P, decreasing = FALSE), ]
  sumstats_top <- sumstats[1:n_snps, ]
  write.table(t(sumstats_top$ID), file = paste0("general_predictors/predictors_top_cis_pval/", protein, ".variants.txt"),
              quote = FALSE, sep = "\n", eol = "\n", row.names = FALSE, col.names = FALSE)
}

