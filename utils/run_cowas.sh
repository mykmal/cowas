#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --partition=agsmall,aglarge,ag2tb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out

# A file that lists pairs of proteins for which to perform the COWAS association test.
# This should be a text file with one tab-separated pair of protein names per line.
PAIRS=pairs/all_protein_pairs.tsv

# File name of the GWAS summary dataset to use
GWAS=data_cleaned/Bellenguez_2022_AD_GWAS.tsv

# Path to a folder containing pre-trained COWAS model weights, as saved by train_cowas.R
WEIGHTS=cowas_weights

# Base name of the genotype data to use as an LD reference panel (in PLINK 2.0 format)
GENOTYPES=data_cleaned/genotypes_subset_for_AD

# Folder with files listing the variants used as predictors for each protein.
# This is only for filtering the genotypes before loading them into cowas.R, in order to decrease runtime.
PREDICTORS=predictors_sis

# File name to which COWAS test results will be written
OUT_FILE=cowas_results.tsv

# Number of cores to use for parallelization
CORES=32


printf "Runtime parameters:\n"
printf "PAIRS = ${PAIRS}\n"
printf "GWAS = ${GWAS}\n"
printf "WEIGHTS = ${WEIGHTS}\n"
printf "GENOTYPES = ${GENOTYPES}\n"
printf "PREDICTORS = ${PREDICTORS}\n"
printf "OUT_FILE = ${OUT_FILE}\n"
printf "CORES = ${CORES}\n\n"

module load R/4.3.0-openblas

# This variable sets the number of cores for R to use
export OMP_NUM_THREADS=${CORES}

if ( [ ! -f ${OUT_FILE} ] ); then
printf "ID_A\tID_B\tN_REFERENCE\tN_GWAS\t\
THETA_DIRECT_A\tVAR_THETA_DIRECT_A\tPVAL_THETA_DIRECT_A\t\
THETA_DIRECT_B\tVAR_THETA_DIRECT_B\tPVAL_THETA_DIRECT_B\t\
THETA_DIRECT_CO\tVAR_THETA_DIRECT_CO\tPVAL_THETA_DIRECT_CO\t\
THETA_FULL_A\tVAR_THETA_FULL_A\tPVAL_THETA_FULL_A\t\
THETA_FULL_B\tVAR_THETA_FULL_B\tPVAL_THETA_FULL_B\t\
THETA_FULL_CO\tVAR_THETA_FULL_CO\tPVAL_THETA_FULL_CO\t\
FSTAT_FULL\tPVAL_FSTAT_FULL\n" > ${OUT_FILE}
fi

# Loop through all protein pairs in the $PAIRS file
while read -r PROTEIN_A PROTEIN_B ETC; do

if ( [ ! -f ${WEIGHTS}/${PROTEIN_A}-${PROTEIN_B}.weights.rds ] ); then
printf "WARNING: no model weights found for ${PROTEIN_A} and ${PROTEIN_B}. Skipping this pair.\n"
continue
fi

if ( [ ! -f ${PREDICTORS}/${PROTEIN_A}.variants.txt ] || [ ! -f ${PREDICTORS}/${PROTEIN_B}.variants.txt ] ); then
printf "WARNING: list of predictors not found for ${PROTEIN_A} or ${PROTEIN_B}. Skipping this pair.\n"
continue
fi

# Create folder for storing temporary files
COWAS_TEMP=TEMP-${PROTEIN_A}-${PROTEIN_B}
mkdir ${COWAS_TEMP}

# Extract genotypes for variants present in the pre-trained models, to use as an LD reference panel
plink2 --pfile ${GENOTYPES} \
       --silent \
       --threads ${CORES} \
       --extract ${PREDICTORS}/${PROTEIN_A}.variants.txt ${PREDICTORS}/${PROTEIN_B}.variants.txt \
       --export A \
       --export-allele data_cleaned/ukb_alt_alleles.tsv \
       --make-just-pvar cols=maybecm \
       --out ${COWAS_TEMP}/genotypes

if ( [ ! -f ${COWAS_TEMP}/genotypes.raw ] ); then
printf "WARNING: Unable to extract genotype data for ${PROTEIN_A} and ${PROTEIN_B}. Skipping this pair.\n"
rm -rf ${COWAS_TEMP}
continue
fi

# Remove unnecessary columns from the raw PLINK genotype file
cut -f 7- ${COWAS_TEMP}/genotypes.raw > ${COWAS_TEMP}/ld_reference.gmatrix

# And this is where the magic happens!
./cowas.R --protein_a ${PROTEIN_A} \
          --protein_b ${PROTEIN_B} \
          --gwas ${GWAS} \
          --weights ${WEIGHTS} \
          --alleles ${COWAS_TEMP}/genotypes.pvar \
          --ld_reference ${COWAS_TEMP}/ld_reference.gmatrix \
          --out ${OUT_FILE} \
          --cores ${CORES}

rm -rf ${COWAS_TEMP}

done < ${PAIRS}

printf "Done!\n"

