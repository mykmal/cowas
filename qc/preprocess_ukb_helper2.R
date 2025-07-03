library(data.table)

# Load data
covariates_nopc <- fread(file = "data_cleaned/covariates_nopc.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
principal_components <- fread(file = "TEMP_top20_pcs.eigenvec", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)

# Rename the #IID column
setnames(principal_components, "#IID", "IID")

# Merge the previously-formatted covariates with the top 20 PCs
covariates_all <- merge(covariates_nopc, principal_components, by = "IID", sort = FALSE)

# Make sure the rows are sorted
setorder(covariates_all, IID)

# Save the final covariate file
fwrite(covariates_all, file = "data_cleaned/covariates.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

