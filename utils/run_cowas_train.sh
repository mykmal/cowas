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

# A file that lists pairs of proteins for which to train COWAS models.
# This should be a text file with one tab-separated pair of protein names per line.
PAIRS=pairs/all_protein_pairs.tsv

# Folder with files listing the variants to use as predictors for each protein.
# Files should be named <PROTEIN_NAME>.variants.txt and contain one column of variant IDs.
PREDICTORS=predictors_sis

# Base name of the genotype data (in PLINK 2.0 format)
GENOTYPES=data_cleaned/genotypes

# File name of the expression data (in plain text, long format)
EXPRESSION=data_cleaned/proteins.tsv

# File name of the covariate data (in plain text, long format)
COVARIATES=data_cleaned/covariates.tsv

# Folder for storing COWAS weights
OUT_DIR=cowas_weights

# The type of model to fit.
# Valid options are stepwise, ridge, lasso, and elastic_net.
MODEL=elastic_net

# Number of cores to use for parallelization
CORES=32

# Correlation threshold for expression and co-expression prediction models
COR_THRESHOLD=0.03

# -------------------------------------------------------------------------------------------------

printf "Runtime parameters:\n"
printf "PAIRS = ${PAIRS}\n"
printf "PREDICTORS = ${PREDICTORS}\n"
printf "GENOTYPES = ${GENOTYPES}\n"
printf "EXPRESSION = ${EXPRESSION}\n"
printf "COVARIATES = ${COVARIATES}\n"
printf "OUT_DIR = ${OUT_DIR}\n"
printf "MODEL = ${MODEL}\n"
printf "CORES = ${CORES}\n"
printf "COR_THRESHOLD = ${COR_THRESHOLD}\n\n"

module load R/4.3.3-openblas

if ( [ ! -d ${OUT_DIR} ] ); then
mkdir ${OUT_DIR}
fi

if ( [ ! -f ${OUT_DIR}/performance_metrics.tsv ] ); then
printf "ID_A\tID_B\tSAMPLE_SIZE\t\
NFEATURES_A\tCORRELATION_A\t\
NFEATURES_B\tCORRELATION_B\t\
NFEATURES_CO\tCORRELATION_CO\n" > ${OUT_DIR}/performance_metrics.tsv
fi

# Read names of proteins with available expression measurements
read -r MEASURED_PROTEINS < ${EXPRESSION}

# Loop through all protein pairs in the $PAIRS file
while read -r PROTEIN_A PROTEIN_B ETC; do

if ( [[ ${PROTEIN_A} == ${PROTEIN_B} ]] ); then
continue
fi

if ( [[ ${MEASURED_PROTEINS} != *"${PROTEIN_A}"* ]] || [[ ${MEASURED_PROTEINS} != *"${PROTEIN_B}"* ]] ); then
printf "WARNING: expression data not found for ${PROTEIN_A} or ${PROTEIN_B}. Skipping this pair.\n"
continue
fi

if ( [ ! -f ${PREDICTORS}/${PROTEIN_A}.variants.txt ] || [ ! -f ${PREDICTORS}/${PROTEIN_B}.variants.txt ] ); then
printf "WARNING: list of predictors not found for ${PROTEIN_A} or ${PROTEIN_B}. Skipping this pair.\n"
continue
fi

# Create folder for storing temporary files
COWAS_TEMP=${OUT_DIR}/TEMP-${PROTEIN_A}-${PROTEIN_B}
mkdir ${COWAS_TEMP}

# Extract pre-screened variants for each protein and export their genotypes to a text file with 0..2 coding
plink2 --pfile ${GENOTYPES} \
       --silent \
       --threads ${CORES} \
       --extract ${PREDICTORS}/${PROTEIN_A}.variants.txt \
       --export A \
       --export-allele data_cleaned/ukb_alt_alleles.tsv \
       --out ${COWAS_TEMP}/${PROTEIN_A}
plink2 --pfile ${GENOTYPES} \
       --silent \
       --threads ${CORES} \
       --extract ${PREDICTORS}/${PROTEIN_B}.variants.txt \
       --export A \
       --export-allele data_cleaned/ukb_alt_alleles.tsv \
       --out ${COWAS_TEMP}/${PROTEIN_B}

if ( [ ! -f ${COWAS_TEMP}/${PROTEIN_A}.raw ] || [ ! -f ${COWAS_TEMP}/${PROTEIN_B}.raw ] ); then
printf "WARNING: Unable to extract genotype data for ${PROTEIN_A} or ${PROTEIN_B}. Skipping this pair.\n"
rm -rf ${COWAS_TEMP}
continue
fi

# Remove unnecessary columns from the raw PLINK genotype files
cut -f 2,7- ${COWAS_TEMP}/${PROTEIN_A}.raw > ${COWAS_TEMP}/${PROTEIN_A}.gmatrix
cut -f 2,7- ${COWAS_TEMP}/${PROTEIN_B}.raw > ${COWAS_TEMP}/${PROTEIN_B}.gmatrix

# And this is where the magic happens!
./cowas_train.R --protein_a ${PROTEIN_A} \
                --protein_b ${PROTEIN_B} \
                --genotypes_a ${COWAS_TEMP}/${PROTEIN_A}.gmatrix \
                --genotypes_b ${COWAS_TEMP}/${PROTEIN_B}.gmatrix \
                --expression ${EXPRESSION} \
                --covariates ${COVARIATES} \
                --out_folder ${OUT_DIR} \
                --model ${MODEL} \
                --cores ${CORES} \
                --cor_threshold ${COR_THRESHOLD} \
                --rank_normalize TRUE

rm -rf ${COWAS_TEMP}

done < ${PAIRS}

printf "Done!\n"

