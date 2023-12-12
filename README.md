# Overview of COWAS

This repository provides a set of tools for performing co-expression wide association studies (COWAS). The goal of COWAS is to test for association between disease and the genetically regulated component of gene or protein co-expression. Here we will always refer to protein expression, but all of the concepts and tools equally apply to gene expression as well.

COWAS is run on one pair of proteins at a time. First, three models are trained on an individual-level reference dataset. One model predicts the expression of the first protein, another model predicts the expression of the second protein, and the third model predicts the co-expression of the two proteins. Then the fitted model weights are used to impute the expression level of each protein and their co-expression into summary-level disease GWAS data. Finally, disease status is jointly tested for association with the imputed expression and co-expression levels. By considering the effect of protein-protein interactions on disease in addition to direct effects, COWAS can identify novel disease-relevant proteins and aid in the interpretation of GWAS findings.

# Software documentation

## Installation

1. Download and unpack the COWAS repository from GitHub.
```bash
wget https://github.com/mykmal/cowas/archive/refs/heads/main.zip
unzip cowas-main.zip && rm cowas-main.zip
mv cowas-main cowas && cd cowas
```
2. Launch R and install the packages optparse, data.table, and glmnet. We used R 4.3.0, optparse 1.7.3, data.table 1.14.8, and glmnet 4.1.8.
```R
install.packages(c("optparse", "data.table", "glmnet"))
```
3. Download PLINK 2.00 and place it in a directory on your $PATH. We used PLINK v2.00a6LM AVX2 AMD (23 Nov 2023).
```bash
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_amd_avx2_20231123.zip
unzip plink2_linux_amd_avx2_20231123.zip && rm plink2_linux_amd_avx2_20231123.zip
sudo mv plink2 /usr/local/bin/
```

## Stage 1: Training co-expression imputation models

## Stage 2: Testing for association between co-expression and disease

# Appendix: Applying COWAS to UK Biobank data

This appendix describes how to obtain and prepare the data we used in our paper.

## UK Biobank downloads

To obtain individual-level data from the UK Biobank you will need to submit an application through the [Access Management System (AMS)](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access). After your application is approved, follow the UK Biobank documentation to download, unpack, and extract the following data fields:

* Data-Field 31 within Category 100094 (self-reported sex)
* Data-Field 34 within Category 100094 (year of birth)
* Data-Field 54 within Category 100024 (UK Biobank assessment center)
* Data-Field 22000 within Category 100313 (genotype measurement batch)
* Data-Field 22006 within Category 100313 (indicator for individuals who self-identified as White British and have very similar genetic ancestry)
* Data-Field 22020 within Category 100313 (indicator for unrelated, high-quality samples used in PCA calculation)
* Data-Field 22828 within Category 100319 (WTCHG imputed genotypes)
* Data-Field 30900 within Category 1839 (UKB-PPP proteomics data)

Next, use PLINK to convert the genotype data to PLINK 2 binary format. Our scripts assume that the genotype data is stored in chromosome-specific files with filenames `ukb_chr<CHR>.pgen` + `ukb_chr<CHR>.psam` + `ukb_chr<CHR>.pvar`. Furthermore, we assume that the proteomics data is stored in a tab-separated file named `olink_data.tsv` beginning with a header line and containing the following four columns: eid, ins_index, protein_id, result. Finally, we assume that the remaining data fields listed above are stored in a single tab-separated file named `ukb_main_dataset.tsv` beginning with a header line and containing the following ten columns: f.eid, f.31.0.0, f.34.0.0, f.54.0.0, f.54.1.0, f.54.2.0, f.54.3.0, f.22000.0.0, f.22006.0.0, f.22020.0.0. Copy all of these files to the folder `cowas/raw`.

## Protein annotations

To facilitate reproducibility, we provide an annotation file for the UK Biobank proteomics data (located at `util/olink_annotations.tsv` in this repository). We created the annotation file by first merging UK Biobank [Resource 1013](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1013) with UK Biobank [Data-Coding 143](https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=143) by gene name. We noticed that one entry (ukb_code 1912) was missing a UniProt ID, so we filled in the correct value from UniProt KB. Next, we exported an annotation file from [Ensembl BioMart](https://grch37.ensembl.org/biomart/martview) with the attributes Gene name, Chromosome/scaffold name, Gene start (bp), Gene end (bp), and UniProtKB/Swiss-Prot ID. After removing entries with patched scaffolds from the BioMart file, we merged it with our annotation file by UniProt ID. Two entries (ukb_code 3 and 163) had multiple matches, so we selected the correct one by comparing gene names. Additionally, we manually looked up the GRCh37 genomic coordinates for 39 entries that were not present in the BioMart file using GeneCards. For 13 entries with multiple genes separated by an underscore and located on the same chromosome (e.g. BOLA2_BOLA2B), we took the union of their genomic regions. Two entries (ukb_code 876 and 1364) correspond to two genes each located on separate chromosomes, so for these we used the coordinates of the second gene. Finally, two other entries (ukb_code 1500 and 2157) correspond to genes for which no mapping exists on the GRCh37 assembly, so we only recorded their chromosomes.

## Alzheimer's disease GWAS

We used summary statistics data from the Alzheimer's disease GWAS published by Bellenguez et al. (2022). The summary-level associations can be downloaded from <https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/>. Place the unpacked file `GCST90027158_buildGRCh38.tsv` in `cowas/raw`.

## Data processing script

Run the shell script `util/preprocess_ukb.sh` from within the main COWAS folder to process the downloaded data. This script will perform the following data wrangling and quality control steps:

1. Create the file `phenotypes/protein_pairs.tsv` listing all possible pairs of proteins that have autosomal genomic coordinates available, with one pair per row.
2. Subset the main dataset to obtain a set of high-quality, unrelated, White British samples. Then subset the genotype data, proteomics data, and covariate data to a common set of samples.
3. Remove variants from the UK Biobank genotype data that are multiallelic, have any missingness, fail a Hardy-Weinberg equilibrium test, lack an rsID, or are palindromic. After this, subset the UK Biobank genotype data and the AD GWAS data to a common set of variants. The filtered genotype and GWAS data are saved to the folders `cowas/genotypes` and `cowas/gwas`, respectively.
4. Compute the top 20 genetic principal components from the quality-controlled genotype data. Following best practices, before computing PCs we remove all indels, SNPs in regions of long-range LD, SNPs with MAF < 0.01, and SNPs with missingness > 0.015. We also prune the remaining SNPs to r^2 < 0.1 with a 1000 bp window and a step size of 80 bp.
5. Create the file `phenotypes/covariates.tsv` with one row per sample and 45 columns for sample ID, sex, year of birth, UKB assessment center (coded as 21 binary dummy variables), genotyping array (binary), and the top 20 genetic PCs. Similarly, create the file `phenotypes/proteins.tsv` with one row per sample and one column per protein.

The `cowas/raw` folder can be deleted after the preprocessing script successfully finishes.

