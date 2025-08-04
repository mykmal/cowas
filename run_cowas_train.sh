#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --partition=msismall,msilarge,msibigmem,msilong
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out

# A file that lists pairs of proteins for which to train COWAS models.
# This should be a text file with one tab-separated pair of protein names per line.
PAIRS=pairs/autosomal_hippie_pairs.tsv

# Folder with files listing the variants to use as predictors for each protein.
# Files should be named <PROTEIN_NAME>.variants.txt and contain a single column of variant IDs.
PREDICTORS=gwas_specific_predictors/predictors_top_cis_beta_pd

# Base name of the genotype data (in PLINK 2.0 format)
GENOTYPES=data_cleaned/genotypes_subset_for_PD

# A file that specifies the counted allele (i.e., the ALT allele) for each variant.
# This should be a text file with one variant per row and two tab-separated columns
# providing the variant ID and the allele that should be counted.
ALT_ALLELES=data_cleaned/ukb_alt_alleles.tsv

# File name of the expression data (in plain text, with samples in rows and proteins in columns)
EXPRESSION=data_cleaned/proteins.tsv

# File name of the covariate data (in plain text, with samples in rows and covariates in columns)
COVARIATES=data_cleaned/covariates.tsv

# Folder for storing COWAS weights
OUT_FOLDER=gwas_specific_weights/pd_weights/weights_cis_beta_elnet_product_pd

# The type of model to fit.
# Valid options are stepwise, ridge, lasso, and elastic_net.
MODEL=elastic_net

# Number of cores to use for parallelization
CORES=32

# Correlation threshold for expression and co-expression imputation models
COR_THRESHOLD=0.03

# Should expression levels be normalized using a rank-based inverse normal transformation?
# If set to FALSE, then expression levels will simply be centered and scaled.
RANK_NORMALIZE=TRUE

# Should product-based COWAS models be trained? If set to FALSE, then residual-based COWAS
# models will be trained instead. In product-based COWAS, the outcome of the co-expression
# model is the product of observed expression levels. In residual-based COWAS, the outcome
# of the co-expression model is the product of single-protein model residuals. See our
# paper for details.
PRODUCT_BASED=TRUE

# -------------------------------------------------------------------------------------------------

printf "Runtime parameters:\n"
printf "PAIRS = ${PAIRS}\n"
printf "PREDICTORS = ${PREDICTORS}\n"
printf "GENOTYPES = ${GENOTYPES}\n"
printf "ALT_ALLELES = ${ALT_ALLELES}\n"
printf "EXPRESSION = ${EXPRESSION}\n"
printf "COVARIATES = ${COVARIATES}\n"
printf "OUT_FOLDER = ${OUT_FOLDER}\n"
printf "MODEL = ${MODEL}\n"
printf "CORES = ${CORES}\n"
printf "COR_THRESHOLD = ${COR_THRESHOLD}\n"
printf "RANK_NORMALIZE = ${RANK_NORMALIZE}\n"
printf "PRODUCT_BASED = ${PRODUCT_BASED}\n\n"

module load R/4.4.2-openblas-rocky8

if ( [ ! -d ${OUT_FOLDER} ] ); then
mkdir -p ${OUT_FOLDER}
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

# Create a folder for storing temporary files
COWAS_TEMP=${OUT_FOLDER}/TEMP-${PROTEIN_A}-${PROTEIN_B}
mkdir ${COWAS_TEMP}

# Extract pre-screened variants for each protein and export their ALT allele dosages to a text file with 0..2 coding
plink2 --pfile ${GENOTYPES} \
       --silent \
       --threads ${CORES} \
       --extract ${PREDICTORS}/${PROTEIN_A}.variants.txt \
       --export A \
       --export-allele ${ALT_ALLELES} \
       --out ${COWAS_TEMP}/${PROTEIN_A}
plink2 --pfile ${GENOTYPES} \
       --silent \
       --threads ${CORES} \
       --extract ${PREDICTORS}/${PROTEIN_B}.variants.txt \
       --export A \
       --export-allele ${ALT_ALLELES} \
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
                --out_folder ${OUT_FOLDER} \
                --model ${MODEL} \
                --cores ${CORES} \
                --cor_threshold ${COR_THRESHOLD} \
                --rank_normalize ${RANK_NORMALIZE} \
                --product_based ${PRODUCT_BASED}

rm -rf ${COWAS_TEMP}

done < ${PAIRS}

printf "Done!\n"

