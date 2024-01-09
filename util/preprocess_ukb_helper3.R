library(data.table)

covariates_nopc <- fread(file = "phenotypes/covariates_nopc.tsv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
principal_components <- fread(file = "TEMP_top20_pcs.eigenvec", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)

covariates_all <- merge(covariates_nopc, principal_components, by = "IID", sort = FALSE)

fwrite(covariates_all, file = "phenotypes/covariates.tsv", quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

