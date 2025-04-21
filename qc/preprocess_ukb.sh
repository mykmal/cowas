#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=124gb
#SBATCH --time=96:00:00
#SBATCH --partition=msismall,msilarge,msibigmem,msilong
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out

module load R/4.4.2-openblas-rocky8

mkdir pairs
mkdir data_cleaned

printf "ukb_chr1\nukb_chr2\nukb_chr3\nukb_chr4\nukb_chr5\n\
ukb_chr6\nukb_chr7\nukb_chr8\nukb_chr9\nukb_chr10\n\
ukb_chr11\nukb_chr12\nukb_chr13\nukb_chr14\nukb_chr15\n\
ukb_chr16\nukb_chr17\nukb_chr18\nukb_chr19\nukb_chr20\n\
ukb_chr21\nukb_chr22\n" > TEMP_plink_files.txt

plink2 --pmerge-list TEMP_plink_files.txt pfile \
       --threads 30 \
       --memory 120000 \
       --pmerge-list-dir data_raw \
       --out TEMP_merged

plink2 --pfile TEMP_merged \
       --threads 30 \
       --memory 120000 \
       --mind 0.01 \
       --make-pgen \
       --out TEMP_1

Rscript --vanilla qc/preprocess_ukb_helper1.R

plink2 --pfile TEMP_1 \
       --threads 30 \
       --memory 120000 \
       --keep TEMP_high_quality_samples.txt \
       --make-pgen \
       --out TEMP_2

plink2 --pfile TEMP_2 \
       --threads 30 \
       --memory 120000 \
       --geno 0.1 \
       --make-pgen \
       --out TEMP_3

plink2 --pfile TEMP_3 \
       --threads 30 \
       --memory 120000 \
       --min-ac 100 \
       --make-pgen \
       --out TEMP_4

plink2 --pfile TEMP_4 \
       --threads 30 \
       --memory 120000 \
       --min-af 0.01 \
       --make-pgen \
       --out TEMP_5

plink2 --pfile TEMP_5 \
       --threads 30 \
       --memory 120000 \
       --hwe 1e-15 midp \
       --nonfounders \
       --make-pgen \
       --out TEMP_6

awk -v FS='\t' '$3 !~ /^rs/ {print $3}' TEMP_6.pvar > TEMP_no_rsid.txt

plink2 --pfile TEMP_6 \
       --threads 30 \
       --memory 120000 \
       --exclude TEMP_no_rsid.txt \
       --make-pgen \
       --out TEMP_7

awk -v FS='\t' '($4 == "A" && $5 == "T") || ($4 == "T" && $5 == "A") || ($4 == "C" && $5 == "G") || ($4 == "G" && $5 == "C") {print $3}' TEMP_7.pvar > TEMP_palindromic.txt

plink2 --pfile TEMP_7 \
       --threads 30 \
       --memory 120000 \
       --exclude TEMP_palindromic.txt \
       --make-pgen \
       --out TEMP_8

plink2 --pfile TEMP_8 \
       --threads 30 \
       --memory 120000 \
       --indep-pairwise 1000 100 0.8 \
       --out TEMP_indep_variants

plink2 --pfile TEMP_8 \
       --threads 30 \
       --memory 120000 \
       --extract TEMP_indep_variants.prune.in \
       --make-pgen psam-cols=sex \
       --out data_cleaned/genotypes

rm TEMP*

# Extract ALT alleles, in order to keep them consistent across COWAS runs
awk -v FS="\t" -v OFS="\t" '!/^#/ {print $3,$5}' data_cleaned/genotypes.pvar > data_cleaned/ukb_alt_alleles.tsv

# Long-range LD regions in GRCh37 are from table S12 of https://www.biorxiv.org/content/10.1101/166298v1
printf "1 48000000 52000000\n\
2 86000000 100500000\n\
2 134500000 138000000\n\
2 183000000 190000000\n\
3 47500000 50000000\n\
3 83500000 87000000\n\
3 89000000 97500000\n\
5 44000000 51500000\n\
5 98000000 100500000\n\
5 129000000 132000000\n\
5 135500000 138500000\n\
6 25000000 33500000\n\
6 57000000 64000000\n\
6 140000000 142500000\n\
7 55000000 66000000\n\
8 8000000 12000000\n\
8 43000000 50000000\n\
10 37000000 43000000\n\
11 45000000 57000000\n\
11 87500000 90500000\n\
12 33000000 40000000\n\
12 109500000 112000000\n\
20 32000000 34500000\n" > TEMP_long_range_ld.txt

plink2 --pfile data_cleaned/genotypes \
       --threads 30 \
       --memory 120000 \
       --exclude bed1 TEMP_long_range_ld.txt \
       --make-pgen \
       --out TEMP_PCA_1

plink2 --pfile TEMP_PCA_1 \
       --threads 30 \
       --memory 120000 \
       --indep-pairwise 1000 100 0.1 \
       --out TEMP_PCA_indep_variants

plink2 --pfile TEMP_PCA_1 \
       --threads 30 \
       --memory 120000 \
       --extract TEMP_PCA_indep_variants.prune.in \
       --make-pgen \
       --out TEMP_PCA_2

plink2 --pfile TEMP_PCA_2 \
       --threads 30 \
       --memory 120000 \
       --pca 20 \
       --out TEMP_top20_pcs

Rscript --vanilla qc/preprocess_ukb_helper2.R

rm data_cleaned/covariates_nopc.tsv

Rscript --vanilla qc/preprocess_ukb_helper3.R

plink2 --pfile data_cleaned/genotypes \
       --threads 30 \
       --memory 120000 \
       --extract TEMP_mutual_LDL_variants.txt \
       --make-pgen psam-cols=sex \
       --out data_cleaned/genotypes_subset_for_LDL

plink2 --pfile data_cleaned/genotypes \
       --threads 30 \
       --memory 120000 \
       --extract TEMP_mutual_AD_variants.txt \
       --make-pgen psam-cols=sex \
       --out data_cleaned/genotypes_subset_for_AD

plink2 --pfile data_cleaned/genotypes \
       --threads 30 \
       --memory 120000 \
       --extract TEMP_mutual_AD_IGAP_variants.txt \
       --make-pgen psam-cols=sex \
       --out data_cleaned/genotypes_subset_for_AD_IGAP

plink2 --pfile data_cleaned/genotypes \
       --threads 30 \
       --memory 120000 \
       --extract TEMP_mutual_PD_variants.txt \
       --make-pgen psam-cols=sex \
       --out data_cleaned/genotypes_subset_for_PD

rm TEMP*

