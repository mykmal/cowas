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

mkdir pqtls

# Perform association testing
plink2 --pfile data_cleaned/ukb_filtered \
       --no-psam-pheno \
       --threads 128 \
       --memory 420000 \
       --pheno data_cleaned/proteins.tsv \
       --not-pheno GLIPR1 \
       --pheno-quantile-normalize \
       --covar data_cleaned/covariates.tsv \
       --covar-variance-standardize \
       --glm omit-ref hide-covar cols=chrom,pos,beta,p \
       --vif 100000 \
       --out TEMP

# Rename and move the results files
for FILENAME in *.glm.linear; do
    NEWFILE=$(echo ${FILENAME} | sed 's/TEMP.//' | sed 's/.glm.linear/.sumstats.tsv/')
    mv ${FILENAME} pqtls/${NEWFILE}
done

# GLIPR1 only has 106 non-missing samples, and within that subset some covariates are constant or collinear.
# Those covariates must be removed before performing a GWAS for this protein, or else PLINK will crash
plink2 --pfile data_cleaned/ukb_filtered \
       --no-psam-pheno \
       --threads 128 \
       --memory 420000 \
       --pheno data_cleaned/proteins.tsv \
       --pheno-name GLIPR1 \
       --pheno-quantile-normalize \
       --covar data_cleaned/covariates.tsv \
       --not-covar ASSESSMENT_CENTER_11001 ASSESSMENT_CENTER_11002 ASSESSMENT_CENTER_11005 ASSESSMENT_CENTER_11018 ASSESSMENT_CENTER_11020 ASSESSMENT_CENTER_11022 ASSESSMENT_CENTER_11023 \
       --covar-variance-standardize \
       --glm omit-ref hide-covar cols=chrom,pos,beta,p \
       --vif 100000 \
       --out TEMP

mv TEMP.GLIPR1.glm.linear pqtls/GLIPR1.sumstats.tsv

rm TEMP*

