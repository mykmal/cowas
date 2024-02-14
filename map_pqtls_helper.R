library(data.table)

protein <- commandArgs(trailingOnly = TRUE)

# Load expression data for the given protein
expression <- fread(file = "data_cleaned/proteins.tsv", header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE,
                    select = c("IID", protein))

# Drop NAs and make sure the data is numeric
expression <- na.omit(expression)
expression[, (protein) := as.numeric(protein)]

# This offset corresponds to the Blom transform
offset = 0.375

# Get needed quantities
n <- nrow(expression)
ranks <- rank(expression[[protein]], ties.method = "average")

# Perform an inverse-rank normal transformation
expression[, (protein) := stats::qnorm((ranks - offset) / (n - 2 * offset + 1))]

# Save data
fwrite(expression, file = paste("TEMP", protein, "normalized.txt", sep = "_"), quote = FALSE, na = "NA", sep = "\t", eol = "\n", compress = "none")

