library(data.table)

### EXTRACT PREDICTORS IN CIS FROM VARIANTS PRESENT IN A GWAS

sumstat_files <- list.files("pqtl_associations/", pattern = "*.sumstats.tsv")
protein_names <- gsub(".sumstats.tsv", "", sumstat_files, fixed = TRUE)

# Change this line to specify a file containing your GWAS variants
variants <- fread(file = "data_cleaned/genotypes_subset_for_AD.pvar", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, select = "ID", verbose = FALSE)

annotations <- fread(file = "olink_annotations_lifted.tsv", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, select = c("Assay", "chr_hg19", "start_hg19", "end_hg19"), verbose = FALSE)
annotations[, start_hg19 := start_hg19 - 500000]
annotations[, end_hg19 := end_hg19 + 500000]

# Select 100 variants with the highest absolute value betas
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
  
  sumstats <- sumstats[ID %in% variants$ID, ]
  sumstats <- sumstats[(CHROM == chr) & (POS >= start) & (POS <= end), ]
  
  sumstats[, BETA := abs(BETA)]
  sumstats <- sumstats[order(sumstats$BETA, decreasing = TRUE), ]
  if (nrow(sumstats) >= 100) {
    sumstats_top <- sumstats[1:100, ]
  } else {
    message(nrow(sumstats), " ", protein)
    sumstats_top <- sumstats[1:nrow(sumstats), ]
  }
  
  # Specify an informative file name here
  write.table(t(sumstats_top$ID), file = paste0("predictors_top_cis_beta_ad/", protein, ".variants.txt"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
}

