#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=124gb
#SBATCH --time=96:00:00
#SBATCH --partition=agsmall,aglarge,ag2tb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out

# Specify the protein for which to perform association testing
# This variable must match the column name for the protein's expression levels
PROTEIN=protein_1

module load R/4.3.0-openblas

if ( [ ! -d pqtls ] ); then
mkdir pqtls
fi

# Perform association testing
plink2 --pfile data_cleaned/ukb_filtered \
       --no-psam-pheno \
       --threads 30 \
	   --memory 120000 \
	   --pheno data_cleaned/proteins.tsv \
	   --pheno-name ${PROTEIN} \
	   --pheno-quantile-normalize \
	   --covar data_cleaned/covariates.tsv \
	   --covar-variance-standardize \
       --glm omit-ref hide-covar cols=p \
	   --vif 100000 \
	   --pfilter 0.01 \
	   --out TEMP_${PROTEIN}

# Rename and move the results file
mv TEMP_${PROTEIN}.${PROTEIN}.glm.linear pqtls/${PROTEIN}_sumstats.tsv
rm TEMP_${PROTEIN}*

