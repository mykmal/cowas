library(data.table)

# Load data
covariates_nopc <- fread(file = "data_cleaned/covariates_nopc.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
principal_components <- fread(file = "TEMP_top20_pcs.eigenvec", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)

# Remove the FID column from the computed PCs, which we don't need
setnames(principal_components, "#FID", "FID")
principal_components[, FID := NULL]

# Merge the previously-formatted covariates with the top 20 PCs
covariates_all <- merge(covariates_nopc, principal_components, by = "IID", sort = FALSE)

# Save the final covariate file
fwrite(covariates_all, file = "data_cleaned/covariates.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

