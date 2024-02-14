#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=420gb
#SBATCH --time=96:00:00
#SBATCH --partition=agsmall,aglarge,ag2tb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out

module load R/4.3.0-openblas

mkdir pqtls

# Perform association testing
plink2 --pfile data_cleaned/ukb_filtered \
	   --no-psam-pheno \
	   --threads 128 \
	   --memory 418000 \
	   --pheno data_cleaned/proteins.tsv \
	   --pheno-quantile-normalize \
	   --covar data_cleaned/covariates.tsv \
	   --covar-variance-standardize \
	   --glm omit-ref hide-covar cols=p \
	   --vif 100000 \
	   --pfilter 0.01 \
	   --out TEMP

# Rename and move the results files
for NUMBER in {1..2923}; do
mv TEMP.protein_${NUMBER}.glm.linear pqtls/protein_${NUMBER}_sumstats.tsv
done

rm TEMP*

