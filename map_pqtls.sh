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

if ( [ ! -d pqtls ] ); then
mkdir pqtls
fi

# Inverse-rank normalize the protein NPX levels
Rscript --vanilla map_pqtls_helper.R ${PROTEIN}

# Perform association testing
plink2 --pfile data_cleaned/ukb_filtered \
       --threads 30 \
	   --memory 120000 \
	   --pheno TEMP_${PROTEIN}_normalized.txt \
	   --no-psam-pheno \
	   --covar data_cleaned/covariates.tsv \
       --glm omit-ref hide-covar cols=p \
	   --pfilter 0.01 \
	   --out TEMP_${PROTEIN}

# Rename and move the results file
mv TEMP_${PROTEIN}.${PROTEIN}.glm.linear pqtls/${PROTEIN}_sumstats.tsv
rm TEMP_${PROTEIN}*

