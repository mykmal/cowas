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

module load R/4.3.0-openblas

mkdir genotypes
mkdir gwas
mkdir phenotypes

Rscript --vanilla util/preprocess_ukb_helper1.R

for CHROM in {1..22}
do
plink2 --pfile raw/ukb_chr${CHROM} \
       --threads 30 \
       --memory 120000 \
       --keep TEMP_high_quality_samples.txt \
       --make-pgen \
       --out TEMP_LOOP_1

plink2 --pfile TEMP_LOOP_1 \
       --threads 30 \
       --memory 120000 \
       --geno 0 \
       --make-pgen \
       --out TEMP_LOOP_2

plink2 --pfile TEMP_LOOP_2 \
       --threads 30 \
       --memory 120000 \
       --hwe 1e-6 midp \
       --nonfounders \
       --make-pgen \
       --out TEMP_LOOP_3

awk -v FS='\t' '$3 !~ /^rs/ {print $3}' TEMP_LOOP_3.pvar > TEMP_LOOP_no_rsid.txt

plink2 --pfile TEMP_LOOP_3 \
       --threads 30 \
       --memory 120000 \
       --exclude TEMP_LOOP_no_rsid.txt \
       --make-pgen \
       --out TEMP_LOOP_4

awk -v FS='\t' '($4 == "A" && $5 == "T") || ($4 == "T" && $5 == "A") || ($4 == "C" && $5 == "G") || ($4 == "G" && $5 == "C") {print $3}' \
          TEMP_LOOP_4.pvar > TEMP_LOOP_palindromic.txt

plink2 --pfile TEMP_LOOP_4 \
       --threads 30 \
       --memory 120000 \
       --exclude TEMP_LOOP_palindromic.txt \
       --make-pgen \
       --out TEMP_FINAL_CHR${CHROM}

rm TEMP_LOOP*
done

cat TEMP_FINAL_CHR1.pvar TEMP_FINAL_CHR2.pvar TEMP_FINAL_CHR3.pvar TEMP_FINAL_CHR4.pvar TEMP_FINAL_CHR5.pvar \
    TEMP_FINAL_CHR6.pvar TEMP_FINAL_CHR7.pvar TEMP_FINAL_CHR8.pvar TEMP_FINAL_CHR9.pvar TEMP_FINAL_CHR10.pvar \
    TEMP_FINAL_CHR11.pvar TEMP_FINAL_CHR12.pvar TEMP_FINAL_CHR13.pvar TEMP_FINAL_CHR14.pvar TEMP_FINAL_CHR15.pvar \
    TEMP_FINAL_CHR16.pvar TEMP_FINAL_CHR17.pvar TEMP_FINAL_CHR18.pvar TEMP_FINAL_CHR19.pvar TEMP_FINAL_CHR20.pvar \
    TEMP_FINAL_CHR21.pvar TEMP_FINAL_CHR22.pvar > TEMP_all_ukb_variants.pvar

Rscript --vanilla util/preprocess_ukb_helper2.R

for CHROM in {1..22}
do
plink2 --pfile TEMP_FINAL_CHR${CHROM} \
       --threads 30 \
       --memory 120000 \
       --extract TEMP_mutual_variants.txt \
       --make-pgen \
       --out genotypes/ukb_filtered_${CHROM}
done

rm TEMP*

printf "ukb_filtered_1\nukb_filtered_2\nukb_filtered_3\nukb_filtered_4\nukb_filtered_5\n\
ukb_filtered_6\nukb_filtered_7\nukb_filtered_8\nukb_filtered_9\nukb_filtered_10\n\
ukb_filtered_11\nukb_filtered_12\nukb_filtered_13\nukb_filtered_14\nukb_filtered_15\n\
ukb_filtered_16\nukb_filtered_17\nukb_filtered_18\nukb_filtered_19\nukb_filtered_20\n\
ukb_filtered_21\nukb_filtered_22\n" > TEMP_plink_files.txt

plink2 --pmerge-list TEMP_plink_files.txt pfile \
       --threads 30 \
       --memory 120000 \
       --pmerge-list-dir genotypes \
       --out TEMP_merged

plink2 --pfile TEMP_merged \
       --threads 30 \
       --memory 120000 \
       --snps-only just-acgt \
       --make-pgen \
       --out TEMP_STEP1

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

plink2 --pfile TEMP_STEP1 \
       --threads 30 \
       --memory 120000 \
       --exclude bed1 TEMP_long_range_ld.txt \
       --make-pgen \
       --out TEMP_STEP2

plink2 --pfile TEMP_STEP2 \
       --threads 30 \
       --memory 120000 \
       --min-af 0.01 \
       --make-pgen \
       --out TEMP_STEP3

plink2 --pfile TEMP_STEP3 \
       --threads 30 \
       --memory 120000 \
       --indep-pairwise 1000 80 0.1 \
       --out TEMP_indep_variants

plink2 --pfile TEMP_STEP3 \
       --threads 30 \
       --memory 120000 \
       --extract TEMP_indep_variants \
       --make-pgen \
       --out TEMP_STEP4

plink2 --pfile TEMP_STEP4 \
       --threads 30 \
       --memory 120000 \
       --pca 20 scols=sid \
       --out TEMP_top20_pcs

join -t "\t" -j 1 --check-order phenotypes/covariates_nopc.tsv TEMP_top20_pcs.eigenvec > phenotypes/covariates.tsv
rm phenotypes/covariates_nopc.tsv

rm TEMP*

