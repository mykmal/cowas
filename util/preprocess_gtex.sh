#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --time=1:00:00
#SBATCH --tmp=10g
#SBATCH --partition=msismall,msilarge,msibigmem
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out

module load R/4.3.0-openblas

awk -v OFS='\t' 'NR > 6 && $1 ~ /^chr[1-9]/ && $3 == "transcript" {print $16,$10,$1,$4,$5,$14}' \
          raw/gencode.v26.GRCh38.genes.gtf | tr -d ";\"" > annotations.txt

./plink --vcf raw/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
          --const-fid \
          --make-bed \
          --out TEMP1

awk -v FS='\t' -v OFS='\t' '(NR > 11) && ($6 == "3") {print "0",$2}' \
          raw/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt > TEMP_EUR.txt

./plink --bfile TEMP1 \
          --keep TEMP_EUR.txt \
          --make-bed \
          --out TEMP2

./plink --bfile TEMP2 \
          --autosome \
          --make-bed \
          --out TEMP3

./plink --bfile TEMP3 \
          --geno 0 \
          --make-bed \
          --out TEMP4

./plink --bfile TEMP4 \
          --hwe 1e-6 midp \
          --make-bed \
          --out TEMP5

./plink --bfile TEMP5 \
          --maf 0.01 \
          --make-bed \
          --out TEMP6

./plink --bfile TEMP6 \
          --snps-only just-acgt \
          --make-bed \
          --out TEMP7

awk -v FS='\t' -v OFS='\t' 'NR > 1 && length($4) == 1 && length($5) == 1 && $7 != "." {print $1,$7}' \
          raw/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > TEMP_rsids.txt

./plink --bfile TEMP7 \
          --update-name TEMP_rsids.txt \
          --make-bed \
          --out TEMP8

awk -v FS='\t' '$2 !~ /^rs/ {print $2}' TEMP8.bim > TEMP_no_rsid.txt

./plink --bfile TEMP8 \
          --exclude TEMP_no_rsid.txt \
          --make-bed \
          --out TEMP9

awk -v FS='\t' '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C" && $6 == "G") || ($5 == "G" && $6 == "C") {print $2}' \
          TEMP9.bim > TEMP_ambiguous.txt

./plink --bfile TEMP9 \
          --exclude TEMP_ambiguous.txt \
          --make-bed \
          --out TEMP_FINAL

mkdir genotypes

for CHR in {1..22}
do
./plink --bfile TEMP_FINAL \
          --chr ${CHR} \
          --make-bed \
          --out genotypes/dosages_chr${CHR}
done

rm TEMP*

mkdir expression
mkdir covariates
mkdir pairs

Rscript --vanilla util/gtex_scripts.R

