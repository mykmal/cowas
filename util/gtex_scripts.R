tissues <-
  c(
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum",
    "Adrenal_Gland",
    "Artery_Aorta",
    "Artery_Coronary",
    "Artery_Tibial",
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
    "Breast_Mammary_Tissue",
    "Cells_Cultured_fibroblasts",
    "Cells_EBV-transformed_lymphocytes",
    "Colon_Sigmoid",
    "Colon_Transverse",
    "Esophagus_Gastroesophageal_Junction",
    "Esophagus_Mucosa",
    "Esophagus_Muscularis",
    "Heart_Atrial_Appendage",
    "Heart_Left_Ventricle",
    "Kidney_Cortex",
    "Liver",
    "Lung",
    "Minor_Salivary_Gland",
    "Muscle_Skeletal",
    "Nerve_Tibial",
    "Ovary",
    "Pancreas",
    "Pituitary",
    "Prostate",
    "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg",
    "Small_Intestine_Terminal_Ileum",
    "Spleen",
    "Stomach",
    "Testis",
    "Thyroid",
    "Uterus",
    "Vagina",
    "Whole_Blood"
  )

for (tissue in tissues) {
  
  expression <- read.table(file = paste0("raw/expression_matrices/", tissue, ".v8.EUR.normalized_expression.bed.gz"),
                           comment.char = "", check.names = FALSE)
  expression <- expression[-c(1,2,3)]
  expression <- as.data.frame(t(expression))
  expression[1, 1] <- "IID"
  write.table(expression, file = paste0("expression/", tissue, ".expression.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  FID <- matrix(0L, nrow = nrow(expression) - 1, ncol = 1)
  individuals <- cbind(FID, expression[, 1][-1])
  write.table(individuals, file = paste0("expression/", tissue, ".individuals.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  covariates <- read.table(file = paste0("raw/expression_covariates/", tissue, ".v8.EUR.covariates.txt"))
  covariates <- as.data.frame(t(covariates))
  covariates[1, 1] <- "IID"
  write.table(covariates, file = paste0("covariates/", tissue, ".covariates.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  annotations <- read.table(file = "annotations.txt", header = FALSE)
  measured_genes <- expression[1, -1]
  
  for (chr in 1:22) {
    current_annotations <- subset(annotations, V3 == paste0("chr", chr))
    current_genes <- intersect(current_annotations$V2, measured_genes)
    pairs <- t(combn(current_genes, 2))
    write.table(pairs, file = paste0("pairs/", tissue, "_chr", chr, ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}

