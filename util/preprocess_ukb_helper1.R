library(data.table)

# Load data
ukb_main_dataset <- fread(file = "raw/ukb_main_dataset.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
olink_data <- fread(file = "raw/olink_data.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
genotyped_samples <- fread(file = "raw/ukb_chr1.psam", header = TRUE, select = 2)
olink_annotations <- fread(file = "util/olink_annotations.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)

# Create a file listing pairs of proteins to analyze
autosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")
olink_annotations <- olink_annotations[chr_hg19 %in% autosomes]
olink_annotations <- na.omit(olink_annotations)
protein_pairs <- expand.grid(olink_annotations$ukb_code, olink_annotations$ukb_code)
write.table(protein_pairs, file = "phenotypes/protein_pairs.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Make sure the columns we will perform operations on are of type integer
ukb_main_dataset[, f.54.0.0 := as.integer(f.54.0.0)]
ukb_main_dataset[, f.22000.0.0 := as.integer(f.22000.0.0)]
ukb_main_dataset[, f.22006.0.0 := as.integer(f.22006.0.0)]
ukb_main_dataset[, f.22020.0.0 := as.integer(f.22020.0.0)]

# Subset to high-quality, unrelated, White British samples
ukb_main_dataset <- ukb_main_dataset[f.22020.0.0 == 1 & f.22006.0.0 == 1]

# Reformat the main dataset
ukb_main_dataset <- ukb_main_dataset[, c("f.54.1.0", "f.54.2.0", "f.54.3.0", "f.22006.0.0", "f.22020.0.0") := NULL]
ukb_main_dataset <- na.omit(ukb_main_dataset)
setnames(ukb_main_dataset, c("IID", "SEX", "BIRTH_YEAR", "ASSESSMENT_CENTER", "MEASUREMENT_BATCH"))

# Create dummy variables for ASSESSMENT_CENTER and MEASUREMENT_BATCH
center_ids <- unique(ukb_main_dataset$ASSESSMENT_CENTER)[-1]
ukb_main_dataset[, (paste0("ASSESSMENT_CENTER_", center_ids)) := lapply(center_ids, function(x) ifelse(ASSESSMENT_CENTER == x, 1, 0))]
ukb_main_dataset[, AXIOM_ARRAY := ifelse(MEASUREMENT_BATCH > 0, 1, 0)]
ukb_main_dataset[, c("ASSESSMENT_CENTER", "MEASUREMENT_BATCH") := NULL]

# Reformat the proteomics data
olink_data[, ins_index := NULL]
olink_data <- na.omit(olink_data)
olink_data <- dcast(olink_data, eid ~ protein_id, value.var = "result")
setnames(olink_data, gsub("result_", "", names(olink_data)))
setnames(olink_data, "eid", "IID")

# Subset all datasets to a common set of samples
common_samples <- Reduce(intersect,
                         list(ukb_main_dataset$IID,
                              olink_data$IID,
                              genotyped_samples$IID))
ukb_main_dataset <- ukb_main_dataset[common_samples]
olink_data <- olink_data[common_samples]

# Put the samples in a format compatible with PLINK --keep
plink_keep_samples <- as.data.frame(cbind(common_samples, common_samples))
names(plink_keep_samples) <- c("#FID", "IID")

# Save the formatted datasets
fwrite(ukb_main_dataset, file = "phenotypes/covariates_nopc.tsv", quote = FALSE, sep = "\t", eol = "\n", compress = "none")
fwrite(olink_data, file = "phenotypes/proteins.tsv", quote = FALSE, sep = "\t", eol = "\n", compress = "none")
write.table(plink_keep_samples, file = "TEMP_high_quality_samples.txt", quote = FALSE, sep = "\t", row.names = FALSE)

