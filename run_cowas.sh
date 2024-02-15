#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=96:00:00
#SBATCH --partition=agsmall,aglarge,ag2tb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out

# A file that lists pairs of genes or proteins for which to perform COWAS
# This should be a text file with one tab-separated pair of gene/protein names per line
PAIRS=pairs/all_protein_pairs.tsv

# File name for COWAS results
OUT_FILE=all_proteins.tsv

# Directory for storing the COWAS output file
OUT_DIR=output_stepwise

# The type of model to fit
# Valid options are stepwise, ridge, and elastic_net
MODEL=stepwise

# R^2 threshold for expression and co-expression prediction models
R2_THRESHOLD=0.001

# Number of cross-validation folds to use for training expression imputation models
# and assessing their predictive performance
CV_FOLDS=10


printf "Runtime parameters:\n"
printf "PAIRS = ${PAIRS}\n"
printf "OUT_FILE = ${OUT_FILE}\n"
printf "OUT_DIR = ${OUT_DIR}\n"
printf "MODEL = ${MODEL}\n"
printf "R2_THRESHOLD = ${R2_THRESHOLD}\n"
printf "CV_FOLDS = ${CV_FOLDS}\n\n"

module load R/4.3.0-openblas

if ( [ ! -d ${OUT_DIR} ] ); then
mkdir ${OUT_DIR}
fi

if [ ! -d ${OUT_DIR}/temp ]; then
mkdir ${OUT_DIR}/temp
fi

if ( [ ! -f ${OUT_DIR}/${OUT_FILE} ] ); then
printf "ID_A\tID_B\tN_REFERENCE\tN_GWAS\t\
NFEATURES_A\tR2PRED_A\t\
NFEATURES_B\tR2PRED_B\t\
NFEATURES_CO\tR2PRED_CO\t\
THETA_DIRECT_A\tVAR_THETA_DIRECT_A\tPVAL_THETA_DIRECT_A\t\
THETA_DIRECT_B\tVAR_THETA_DIRECT_B\tPVAL_THETA_DIRECT_B\t\
THETA_DIRECT_CO\tVAR_THETA_DIRECT_CO\tPVAL_THETA_DIRECT_CO\t\
THETA_FULL_A\tVAR_THETA_FULL_A\tPVAL_THETA_FULL_A\t\
THETA_FULL_B\tVAR_THETA_FULL_B\tPVAL_THETA_FULL_B\t\
THETA_FULL_CO\tVAR_THETA_FULL_CO\tPVAL_THETA_FULL_CO\t\
FSTAT_FULL\tPVAL_FSTAT_FULL\n" > ${OUT_DIR}/${OUT_FILE}
fi

# Read names of proteins with available expression measurements
read -r MEASURED_PROTEINS < data_cleaned/proteins.tsv

# Loop through all protein pairs in the $PAIRS file
while read -r PROTEIN_A PROTEIN_B ETC; do

if ( [[ ${PROTEIN_A} == ${PROTEIN_B} ]] ); then
continue
fi

if ( [[ ${MEASURED_PROTEINS} != *"${PROTEIN_A}"* ]] || [[ ${MEASURED_PROTEINS} != *"${PROTEIN_B}"* ]] ); then
printf "WARNING: expression data not found for ${PROTEIN_A} or ${PROTEIN_B}. Skipping this pair.\n"
continue
fi

# Create folder for storing temporary files
COWAS_TEMP=${OUT_DIR}/temp/${PROTEIN_A}_${PROTEIN_B}
mkdir ${COWAS_TEMP}

# Extract pre-screened variants for each protein, then create a pvar file listing them
# and export their genotypes to a text file with 0..2 coding
plink2 --pfile data_cleaned/ukb_filtered \
       --silent \
	   --extract variants/${PROTEIN_A}_features.txt \
	   --export A \
	   --export-allele data_cleaned/ukb_alt_alleles.tsv \
	   --make-just-pvar cols=maybecm \
	   --out ${COWAS_TEMP}/${PROTEIN_A}
plink2 --pfile data_cleaned/ukb_filtered \
       --silent \
	   --extract variants/${PROTEIN_B}_features.txt \
	   --export A \
	   --export-allele data_cleaned/ukb_alt_alleles.tsv \
	   --make-just-pvar cols=maybecm \
	   --out ${COWAS_TEMP}/${PROTEIN_B}

if ( [ ! -f ${COWAS_TEMP}/${PROTEIN_A}.raw ] || [ ! -f ${COWAS_TEMP}/${PROTEIN_B}.raw ] ); then
printf "WARNING: Unable to extract genotype data for ${PROTEIN_A} or ${PROTEIN_B}. Skipping this pair.\n"
rm -rf ${COWAS_TEMP}
continue
fi

# Remove unnecessary columns from the raw PLINK genotype files
cut -f 2,7- ${COWAS_TEMP}/${PROTEIN_A}.raw > ${COWAS_TEMP}/${PROTEIN_A}.gmatrix
cut -f 2,7- ${COWAS_TEMP}/${PROTEIN_B}.raw > ${COWAS_TEMP}/${PROTEIN_B}.gmatrix

# Create a single pvar file with variants for both proteins
# This will duplicate variants that appear in both files, but cowas.R takes care of that later
cat ${COWAS_TEMP}/${PROTEIN_A}.pvar ${COWAS_TEMP}/${PROTEIN_B}.pvar > ${COWAS_TEMP}/both.snps

# And this is where the magic happens!
./cowas.R --protein_a ${PROTEIN_A} \
          --protein_b ${PROTEIN_B} \
		  --genotypes_a ${COWAS_TEMP}/${PROTEIN_A}.gmatrix \
		  --genotypes_b ${COWAS_TEMP}/${PROTEIN_B}.gmatrix \
		  --snps ${COWAS_TEMP}/both.snps \
		  --expression data_cleaned/proteins.tsv \
		  --covariates data_cleaned/covariates.tsv \
		  --gwas data_cleaned/Bellenguez_2022_AD_gwas.tsv \
		  --out ${OUT_DIR}/${OUT_FILE} \
		  --model ${MODEL} \
		  --r2_threshold ${R2_THRESHOLD} \
		  --cv_folds ${CV_FOLDS}

rm -rf ${COWAS_TEMP}

done < ${PAIRS}

printf "Done!\n"

