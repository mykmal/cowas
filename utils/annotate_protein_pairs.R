# This script provides the code we used to identify pairs of proteins in the UK Biobank (UKB)
# that have some prior evidence of interactions and are coded by autosomal genes.

library(tidyr)

# Load the UKB protein annotation file, which we downloaded from https://www.synapse.org/Synapse:syn52364558.
# The protein ANP32C on line 1544 had incorrect positions in the original file, so we manually corrected it.
olink_annotations <- read.delim("protein_annotations_derived/1_olink_protein_map_3k_v1_corrected_typo.tsv")

# Extract the position information for each gene
olink_positions <- olink_annotations[, c("chr", "gene_start", "gene_end")]
olink_positions$chr <- paste0("chr", olink_positions$chr)

# Save the position information
write.table(olink_positions, file = "protein_annotations_derived/2_positions_hg38.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Next, upload this file to the UCSC LiftOver web tool (https://genome.ucsc.edu/cgi-bin/hgLiftOver)
# and convert it from hg38 to hg19. We set the minimum ratio of bases that must overlap to 0.1,
# and unchecked the "Keep original positions in output" and "Allow multiple output regions" boxes.

# Load the file that we just lifted over
lifted_positions <- read.delim("protein_annotations_derived/3_positions_hg19_liftover.bed",
                               header = FALSE)
names(lifted_positions) <- c("chr_hg19", "start_hg19", "end_hg19")
lifted_positions$chr_hg19 <- gsub("chr", "", lifted_positions$chr_hg19, fixed = TRUE)

# Add the lifted positions back into the annotation file
olink_annotations <- subset(olink_annotations, select = -c(chr, gene_start, gene_end))
olink_annotations_lifted <- cbind(olink_annotations, lifted_positions)

# Save the results
write.table(olink_annotations_lifted, file = "protein_annotations_derived/4_olink_annotations_lifted.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Load the HIPPIE database, which we downloaded from https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php
hippie_current <- read.delim("protein_annotations_downloaded/hippie_current.txt",
                             header = FALSE)
names(hippie_current) <- c("UniProt_A", "ID_A", "UniProt_B", "ID_B", "ConfidenceScore", "Evidence")

# Some rows have multiple protein names separated by commas, so we split these into separate rows
hippie_split_rows <- separate_longer_delim(hippie_current, cols = "UniProt_A", delim = ",")
hippie_split_rows <- separate_longer_delim(hippie_split_rows, cols = "UniProt_B", delim = ",")

# Save a version of the HIPPIE database with one protein identifier per row
write.table(hippie_split_rows, file = "protein_annotations_derived/5_hippie_uncollapsed_ids.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Create a vector of all protein identifiers in HIPPIE
hippie_unique_proteins <- unique(c(hippie_split_rows$UniProt_A, hippie_split_rows$UniProt_B))

# Save a file with all protein names
write.table(hippie_unique_proteins, file = "protein_annotations_derived/6_hippie_unique_proteins.txt",
            quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)

# Next, use the UniProt ID Mapping web tool (https://www.uniprot.org/id-mapping)
# to map this file from "UniProt / UniProtKB AC/ID" to "UniProt / Gene Name"

# Load the results from the ID Mapping web tool
uniprot_mapping_from_hippie_to_genes <- read.delim("protein_annotations_derived/7_uniprot_mapping_from_hippie_to_genes.tsv")

# Subset the UniProt mapping file to proteins found in UKB
uniprot_mapping_file_subset_to_ukb <- uniprot_mapping_from_hippie_to_genes[uniprot_mapping_from_hippie_to_genes$To %in% c(olink_annotations_lifted$Assay, olink_annotations_lifted$HGNC.symbol), ]

# Save the results
write.table(uniprot_mapping_file_subset_to_ukb, file = "protein_annotations_derived/8_uniprot_mapping_file_subset_to_ukb.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Map from UniProt IDs to gene names for the second protein in each pair
hippie_subset_to_ukb <- merge(hippie_split_rows, uniprot_mapping_file_subset_to_ukb,
                              by.x = "UniProt_B", by.y = "From")
names(hippie_subset_to_ukb)[7] <- "Gene_B"

# Map from UniProt IDs to gene names for the first protein in each pair
hippie_subset_to_ukb <- merge(hippie_subset_to_ukb, uniprot_mapping_file_subset_to_ukb,
                              by.x = "UniProt_A", by.y = "From")
names(hippie_subset_to_ukb)[8] <- "Gene_A"

# Organize and save the results
hippie_subset_to_ukb <- hippie_subset_to_ukb[, c("UniProt_A", "ID_A", "Gene_A", "UniProt_B", "ID_B", "Gene_B", "ConfidenceScore", "Evidence")]
write.table(hippie_subset_to_ukb, file = "protein_annotations_derived/9_hippie_subset_to_ukb.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Subset the annotations to autosomal genes
autosomal_annotations <- olink_annotations_lifted[olink_annotations_lifted$chr_hg19 %in% 1:22, ]

# Subset the HIPPIE table to pairs of autosomal genes
hippie_subset_to_ukb_autosomal <- hippie_subset_to_ukb[hippie_subset_to_ukb$Gene_A %in% c(autosomal_annotations$Assay, autosomal_annotations$HGNC.symbol) &
                                                         hippie_subset_to_ukb$Gene_B %in% c(autosomal_annotations$Assay, autosomal_annotations$HGNC.symbol), ]

# Save the results
write.table(hippie_subset_to_ukb_autosomal, file = "protein_annotations_derived/10_hippie_subset_to_ukb_autosomal.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Add columns for assay names to the autosomal HIPPIE table
hippie_ukb_assays <- hippie_subset_to_ukb_autosomal
hippie_ukb_assays$Assay_A <- ""
hippie_ukb_assays$Assay_B <- ""

# Get the assay names for each protein. We cannot use a simple merge here because sometimes the
# gene names obtained from the UniProt ID Mapping tool match the UKB "Assay" column and sometimes
# they match the "HGNC.symbol" column.
for (i in 1:nrow(hippie_ukb_assays)) {
  # First protein in the pair
  if (hippie_ukb_assays[i, "Gene_A"] %in% autosomal_annotations$HGNC.symbol) {
    hippie_ukb_assays[i, "Assay_A"] <- autosomal_annotations[autosomal_annotations$HGNC.symbol == hippie_ukb_assays[i, "Gene_A"], "Assay"][1]
  } else {
    hippie_ukb_assays[i, "Assay_A"] <- autosomal_annotations[autosomal_annotations$Assay == hippie_ukb_assays[i, "Gene_A"], "Assay"][1]
  }
  
  # Second protein in the pair
  if (hippie_ukb_assays[i, "Gene_B"] %in% autosomal_annotations$HGNC.symbol) {
    hippie_ukb_assays[i, "Assay_B"] <- autosomal_annotations[autosomal_annotations$HGNC.symbol == hippie_ukb_assays[i, "Gene_B"], "Assay"][1]
  } else {
    hippie_ukb_assays[i, "Assay_B"] <- autosomal_annotations[autosomal_annotations$Assay == hippie_ukb_assays[i, "Gene_B"], "Assay"][1]
  }
}

# Save the results
write.table(hippie_ukb_assays, file = "protein_annotations_derived/11_hippie_subset_to_ukb_autosomal_with_assays.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Drop all columns except the assay names, since that's the only information we need now
hippie_ukb_assays_sorted <- hippie_ukb_assays[, c("Assay_A", "Assay_B")]

# Remove pairs of the same protein
hippie_ukb_assays_sorted <- hippie_ukb_assays_sorted[hippie_ukb_assays_sorted$Assay_A != hippie_ukb_assays_sorted$Assay_B, ]

# Sort the columns within each row alphabetically
for (i in 1:nrow(hippie_ukb_assays_sorted)) {
 sorted <- sort(c(hippie_ukb_assays_sorted[i, 1], hippie_ukb_assays_sorted[i, 2]))
 hippie_ukb_assays_sorted[i, 1] <- sorted[1]
 hippie_ukb_assays_sorted[i, 2] <- sorted[2]
}

# Remove duplicated pairs
hippie_ukb_assays_sorted$pasted_ids <- paste0(hippie_ukb_assays_sorted$Assay_A, "_", hippie_ukb_assays_sorted$Assay_B)
hippie_ukb_assays_sorted <- hippie_ukb_assays_sorted[!duplicated(hippie_ukb_assays_sorted$pasted_ids), ]
hippie_ukb_assays_sorted <- subset(hippie_ukb_assays_sorted, select = -pasted_ids)

# Save the final table of protein pairs that will be used to train COWAS models
write.table(hippie_ukb_assays_sorted, file = "pairs/autosomal_hippie_pairs.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

