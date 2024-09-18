library(data.table)

### EXTRACT PREDICTORS FROM ALL VARIANTS

sumstat_files <- list.files("pqtl_associations/", pattern = "*.sumstats.tsv")
protein_names <- gsub(".sumstats.tsv", "", sumstat_files, fixed = TRUE)

# Select 100 variants with the highest absolute value betas
for (protein in protein_names) {
  sumstats <- fread(file = paste0("pqtl_associations/", protein, ".sumstats.tsv"), header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, select = c("ID", "BETA"), verbose = FALSE)
  sumstats[, BETA := abs(BETA)]
  sumstats <- sumstats[order(sumstats$BETA, decreasing = TRUE), ]
  sumstats_top <- sumstats[1:100, ]
  write.table(t(sumstats_top$ID), file = paste0("predictors_top_global_beta/", protein, ".variants.txt"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
}

# Select 100 variants with the lowest P-values
for (protein in protein_names) {
  sumstats <- fread(file = paste0("pqtl_associations/", protein, ".sumstats.tsv"), header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, select = c("ID", "P"), verbose = FALSE)
  sumstats[, P := as.numeric(P)]
  sumstats <- sumstats[order(sumstats$P, decreasing = FALSE), ]
  sumstats_top <- sumstats[1:100, ]
  write.table(t(sumstats_top$ID), file = paste0("predictors_top_global_pval/", protein, ".variants.txt"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
}


### EXTRACT PREDICTORS IN CIS FROM ALL VARIANTS

sumstat_files <- list.files("pqtl_associations/", pattern = "*.sumstats.tsv")
protein_names <- gsub(".sumstats.tsv", "", sumstat_files, fixed = TRUE)
annotations <- fread(file = "olink_annotations_lifted.tsv", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, select = c("Assay", "chr_hg19", "start_hg19", "end_hg19"), verbose = FALSE)

annotations[, start_hg19 := start_hg19 - 500000]
annotations[, end_hg19 := end_hg19 + 500000]

# Select 100 cis-variants with the highest absolute value betas
for (protein in protein_names) {
  annotations_subset <- annotations[Assay == protein, ]
  
  if (nrow(annotations_subset) > 1) {
    annotations_subset <- annotations_subset[1, ]
  }
  
  chr <- annotations_subset$chr_hg19
  if (chr == "chrX") {
    next
  }
  chr <- as.integer(gsub("chr", "", chr))
  start <- annotations_subset$start_hg19
  end <- annotations_subset$end_hg19
  
  sumstats <- fread(file = paste0("pqtl_associations/", protein, ".sumstats.tsv"), header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, select = c("#CHROM", "POS", "ID", "BETA"), verbose = FALSE)
  setnames(sumstats, "#CHROM", "CHROM")
  sumstats <- sumstats[(CHROM == chr) & (POS >= start) & (POS <= end), ]
 
  sumstats[, BETA := abs(BETA)]
  sumstats <- sumstats[order(sumstats$BETA, decreasing = TRUE), ]
  if (nrow(sumstats) >= 100) {
    sumstats_top <- sumstats[1:100, ]
  } else {
    message(nrow(sumstats), " ", protein)
    sumstats_top <- sumstats[1:nrow(sumstats), ]
  }
  
  write.table(t(sumstats_top$ID), file = paste0("predictors_top_cis_beta/", protein, ".variants.txt"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
}

# Select 100 cis-variants with the lowest P-values
for (protein in protein_names) {
  annotations_subset <- annotations[Assay == protein, ]
  
  if (nrow(annotations_subset) > 1) {
    annotations_subset <- annotations_subset[1, ]
  }
  
  chr <- annotations_subset$chr_hg19
  if (chr == "chrX") {
    next
  }
  chr <- as.integer(gsub("chr", "", chr))
  start <- annotations_subset$start_hg19
  end <- annotations_subset$end_hg19
  
  sumstats <- fread(file = paste0("pqtl_associations/", protein, ".sumstats.tsv"), header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, select = c("#CHROM", "POS", "ID", "P"), verbose = FALSE)
  setnames(sumstats, "#CHROM", "CHROM")
  sumstats <- sumstats[(CHROM == chr) & (POS >= start) & (POS <= end), ]
 
  sumstats[, P := as.numeric(P)]
  sumstats <- sumstats[order(sumstats$P, decreasing = FALSE), ]
  if (nrow(sumstats) >= 100) {
    sumstats_top <- sumstats[1:100, ]
  } else {
    message(nrow(sumstats), " ", protein)
    sumstats_top <- sumstats[1:nrow(sumstats), ]
  }
  
  write.table(t(sumstats_top$ID), file = paste0("predictors_top_cis_pval/", protein, ".variants.txt"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
}

