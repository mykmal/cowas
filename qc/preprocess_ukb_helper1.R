library(data.table)

# Load data
ukb_main_dataset <- fread(file = "data_raw/ukb_main_dataset.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE,
                          select = c("f.eid", "f.31.0.0", "f.54.0.0", "f.21003.0.0", "f.22000.0.0", "f.22006.0.0", "f.22020.0.0"))
olink_data <- fread(file = "data_raw/protein_npx_levels.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
olink_coding <- fread(file = "data_raw/coding143.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
genotyped_samples <- fread(file = "TEMP_merged.psam", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE, select = 2)

# Subset to high-quality, unrelated samples of White British ancestry
ukb_main_dataset <- ukb_main_dataset[f.22020.0.0 == 1 & f.22006.0.0 == 1, ]

# Reformat the main dataset
ukb_main_dataset[, c("f.22006.0.0", "f.22020.0.0") := NULL]
ukb_main_dataset <- na.omit(ukb_main_dataset)
setnames(ukb_main_dataset, c("IID", "SEX", "ASSESSMENT_CENTER", "AGE", "MEASUREMENT_BATCH"))

# Keep samples assessed at the initial visit
olink_data <- olink_data[ins_index == 0, ]
olink_data[, ins_index := NULL]
olink_data <- na.omit(olink_data)

# Replace dummy protein identifiers with their gene names
olink_coding[, c("gene_id", "full_name") := tstrsplit(meaning, ";", fixed = TRUE)]
olink_data[, protein_id_replaced := olink_coding[match(olink_data$protein_id, olink_coding$coding), "gene_id"]]

# Convert the proteomic data from long to wide format
olink_data <- dcast(olink_data, eid ~ protein_id_replaced, value.var = "result")
setnames(olink_data, "eid", "IID")

# Remove GLIPR1 because 99.4% of samples failed quality control for this protein
# See https://www.nature.com/articles/s41586-023-06592-6#Sec20
olink_data[, GLIPR1 := NULL]

# Subset all datasets to a common set of samples
common_samples <- Reduce(intersect,
                         list(ukb_main_dataset$IID,
                              olink_data$IID,
                              genotyped_samples$IID))
ukb_main_dataset <- ukb_main_dataset[IID %in% common_samples, ]
olink_data <- olink_data[IID %in% common_samples, ]

# Create interaction variables for age and sex
ukb_main_dataset[, AGE := as.numeric(AGE)]
ukb_main_dataset[, SEX := as.numeric(SEX)]
ukb_main_dataset[, AGE2 := AGE^2]
ukb_main_dataset[, AGE_SEX := AGE * SEX]
ukb_main_dataset[, AGE2_SEX := AGE2 * SEX]

# Create dummy variables for ASSESSMENT_CENTER and MEASUREMENT_BATCH
center_ids <- sort(as.integer(unique(ukb_main_dataset$ASSESSMENT_CENTER)))[-1]
ukb_main_dataset[, paste0("ASSESSMENT_CENTER_", center_ids) := lapply(center_ids, function(x) ifelse(ASSESSMENT_CENTER == x, 1, 0))]
ukb_main_dataset[, AXIOM_ARRAY := ifelse(MEASUREMENT_BATCH > 0, 1, 0)]
ukb_main_dataset[, c("ASSESSMENT_CENTER", "MEASUREMENT_BATCH") := NULL]
sorted_ukb_columns <- c("IID", "AGE", "AGE2", "SEX", "AGE_SEX", "AGE2_SEX", paste0("ASSESSMENT_CENTER_", center_ids), "AXIOM_ARRAY")
ukb_main_dataset <- ukb_main_dataset[, ..sorted_ukb_columns]

# Remove covariates that have NAs or have become constant after the subsetting
for (column in names(ukb_main_dataset)[-1]) {
  if (anyNA(ukb_main_dataset[[column]]) || var(ukb_main_dataset[[column]], na.rm = TRUE) <= 0) {
    ukb_main_dataset[, (column) := NULL]
  }
}

# Remove proteins that have become fully NA or constant after the subsetting.
# Note that here we allow some NAs to remain because none of the individuals have data for all proteins.
for (column in names(olink_data)[-1]) {
  if (all(is.na(olink_data[[column]])) || var(olink_data[[column]], na.rm = TRUE) <= 0) {
    olink_data[, (column) := NULL]
  }
}

# Create a table listing the pairs of proteins to analyze
protein_pairs <- t(combn(names(olink_data)[-1], 2))

# Sort the rows by sample identifiers
setorder(ukb_main_dataset, IID)
setorder(olink_data, IID)

# Put the samples in a format compatible with PLINK's --keep command
plink_keep_samples <- as.data.frame(cbind(ukb_main_dataset$IID, ukb_main_dataset$IID))
names(plink_keep_samples) <- c("#FID", "IID")

# Save the formatted datasets
fwrite(ukb_main_dataset, file = "data_cleaned/covariates_nopc.tsv",
       quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")
fwrite(olink_data, file = "data_cleaned/proteins.tsv",
       quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")
write.table(protein_pairs, file = "pairs/all_protein_pairs.tsv",
            quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(plink_keep_samples, file = "TEMP_high_quality_samples.txt",
            quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = TRUE)

