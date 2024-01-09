library(data.table)

# Load data
ukb_main_dataset <- fread(file = "raw/ukb_main_dataset.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
olink_data <- fread(file = "raw/olink_data.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
genotyped_samples <- fread(file = "raw/ukb_chr1.psam", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE, select = 2)
olink_annotations <- fread(file = "util/olink_annotations.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)

# Create a file listing pairs of proteins to analyze
autosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")
olink_annotations <- olink_annotations[chr_hg19 %in% autosomes, ]
olink_annotations <- na.omit(olink_annotations)
protein_pairs <- expand.grid(olink_annotations$ukb_code, olink_annotations$ukb_code)
write.table(protein_pairs, file = "phenotypes/protein_pairs.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Subset to high-quality, unrelated, White British samples
ukb_main_dataset <- ukb_main_dataset[f.22020.0.0 == 1 & f.22006.0.0 == 1, ]

# Reformat the main dataset
ukb_main_dataset[, c("f.54.1.0", "f.54.2.0", "f.54.3.0", "f.21003.1.0", "f.21003.2.0", "f.21003.3.0", "f.22006.0.0", "f.22020.0.0") := NULL]
ukb_main_dataset <- na.omit(ukb_main_dataset)
setnames(ukb_main_dataset, c("IID", "SEX", "ASSESSMENT_CENTER", "AGE", "MEASUREMENT_BATCH"))

# Create interaction variables for age and sex
ukb_main_dataset[, AGE := as.integer(AGE)]
ukb_main_dataset[, SEX := as.integer(SEX)]
ukb_main_dataset[, AGE2 := AGE^2]
ukb_main_dataset[, AGE_SEX := AGE * SEX]
ukb_main_dataset[, AGE2_SEX := AGE2 * SEX]

# Create dummy variables for ASSESSMENT_CENTER and MEASUREMENT_BATCH
center_ids <- unique(ukb_main_dataset$ASSESSMENT_CENTER)[-1]
ukb_main_dataset[, paste0("ASSESSMENT_CENTER_", center_ids) := lapply(center_ids, function(x) ifelse(ASSESSMENT_CENTER == x, 1, 0))]
ukb_main_dataset[, AXIOM_ARRAY := ifelse(MEASUREMENT_BATCH > 0, 1, 0)]
ukb_main_dataset[, c("ASSESSMENT_CENTER", "MEASUREMENT_BATCH") := NULL]
ukb_main_dataset[, c("IID", "AGE", "AGE2", "SEX", "AGE_SEX", "AGE2_SEX", paste0("ASSESSMENT_CENTER_", center_ids), "AXIOM_ARRAY")]

# Reformat the proteomics data
olink_data <- olink_data[ins_index == 0, ]
olink_data[, ins_index := NULL]
olink_data <- na.omit(olink_data)
olink_data <- dcast(olink_data, eid ~ protein_id, value.var = "result")
all_proteins <- sort(as.integer(names(olink_data)[-1]))
olink_data[, c("eid", all_proteins)]
setnames(olink_data, paste0("protein_", names(olink_data)))
setnames(olink_data, "protein_eid", "IID")

# Subset all datasets to a common set of samples
common_samples <- Reduce(intersect,
                         list(ukb_main_dataset$IID,
                              olink_data$IID,
                              genotyped_samples$IID))
ukb_main_dataset <- ukb_main_dataset[IID %in% common_samples, ]
olink_data <- olink_data[IID %in% common_samples, ]

# Put the samples in a format compatible with PLINK --keep
plink_keep_samples <- as.data.frame(cbind(common_samples, common_samples))
names(plink_keep_samples) <- c("#FID", "IID")

# Save the formatted datasets
fwrite(ukb_main_dataset, file = "phenotypes/covariates_nopc.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")
fwrite(olink_data, file = "phenotypes/proteins.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")
write.table(plink_keep_samples, file = "TEMP_high_quality_samples.txt", quote = FALSE, sep = "\t", row.names = FALSE)

