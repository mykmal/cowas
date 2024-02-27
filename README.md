# Overview of COWAS

This repository provides a set of tools for performing co-expression wide association studies (COWAS). The goal of COWAS is to test for association between disease and the genetically regulated component of gene or protein co-expression. Here we will always refer to protein expression, but all of the concepts and tools equally apply to gene expression as well.

COWAS is run on one pair of proteins at a time. First, three models are trained on an individual-level reference dataset. One model predicts the expression of the first protein, another model predicts the expression of the second protein, and the third model predicts the co-expression of the two proteins. Then the fitted model weights are used to impute the expression level of each protein and their co-expression into summary-level disease GWAS data. Finally, disease status is jointly tested for association with the imputed expression and co-expression levels. By considering the effect of protein-protein interactions on disease in addition to direct effects, COWAS can identify novel disease-relevant proteins and aid in the interpretation of GWAS findings.

# Installation

1. Download and unpack the COWAS repository from GitHub.
```bash
wget https://github.com/mykmal/cowas/archive/refs/heads/main.zip
unzip main.zip && rm main.zip
mv cowas-main cowas && cd cowas
```
2. Launch R and install the required packages optparse and data.table. If you wish to use ridge, lasso, or elastic net models then also install the package glmnet. If you wish to utilize parallel computation in glmnet then also install the package doMC. We used R 4.3.0, optparse 1.7.4, data.table 1.15.0, glmnet 4.1-8, and doMC 1.3.8.
```R
install.packages(c("optparse", "data.table", "glmnet", "doMC"))
```
3. Download PLINK 2.00 and place it in a directory on your PATH. We used PLINK v2.00a6LM AVX2 AMD (5 Feb 2024).
```bash
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_amd_avx2_20240205.zip
unzip plink2_linux_amd_avx2_20240205.zip && rm plink2_linux_amd_avx2_20240205.zip
sudo mv plink2 /usr/local/bin/
```

# Data preparation and QC

This section describes how to obtain and prepare the data we used in our paper.

## UK Biobank downloads

To obtain individual-level data from the UK Biobank you will need to submit an application through the [Access Management System (AMS)](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access). After your application is approved, follow the UK Biobank documentation to download, unpack, and extract the following data fields:

* Data-Field 31 within Category 100094 (self-reported sex)
* Data-Field 54 within Category 100024 (UK Biobank assessment center)
* Data-Field 21003 within Category 100024 (age when attended assessment center)
* Data-Field 22000 within Category 100313 (genotype measurement batch)
* Data-Field 22006 within Category 100313 (indicator for individuals who self-identified as White British and have very similar genetic ancestry)
* Data-Field 22020 within Category 100313 (indicator for unrelated, high-quality samples used in PCA calculation)
* Data-Field 22828 within Category 100319 (WTCHG imputed genotypes)
* Data-Field 30900 within Category 1839 (UKB-PPP proteomics data)

Next, use PLINK to convert the genotype data to PLINK 2 binary format. Our scripts assume that the genotype data is stored in chromosome-specific files with filenames `ukb_chr<CHR>.pgen` + `ukb_chr<CHR>.psam` + `ukb_chr<CHR>.pvar`. Furthermore, we assume that the proteomics data is stored in a tab-separated file named `olink_data.tsv` beginning with a header line and containing the following four columns: eid, ins_index, protein_id, result. Finally, we assume that the remaining data fields listed above are stored in a single tab-separated file named `ukb_main_dataset.tsv` beginning with a header line and containing the following 13 columns: f.eid, f.31.0.0, f.54.0.0, f.54.1.0, f.54.2.0, f.54.3.0, f.21003.0.0, f.21003.1.0, f.21003.2.0, f.21003.3.0, f.22000.0.0, f.22006.0.0, f.22020.0.0. Copy all of these files to a new subfolder named `data_raw` within your main COWAS folder.

## Alzheimer's disease GWAS

We used summary statistics data from the Alzheimer's disease GWAS published by Bellenguez et al. (2022). The summary-level associations can be downloaded from <https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/>. Place the unpacked file `GCST90027158_buildGRCh38.tsv` in `data_raw`.

## Data processing script

Run the shell script `util/preprocess_ukb.sh` from within the main COWAS folder to process the downloaded data. This script will perform the following data wrangling and quality control steps:

1. Create the file `pairs/all_protein_pairs.tsv` listing all possible pairs of proteins that have data available, with one pair per row.
2. Filter the main dataset to obtain a set of high-quality, unrelated, White British samples with per-sample genotyping rate > 99%. Filter the proteomics data to obtain the set of samples assessed at the initial visit. Then subset the genotype data, proteomics data, and covariate data to a common set of samples.
3. Remove variants from the UK Biobank genotype data that have a missingness rate > 10%, have an MAC < 100, have an MAF < 1%, fail a Hardy-Weinberg equilibrium test ($P < 10^{-15}$), lack an rsID, or are palindromic. Next, prune the remaining variants to R^2 < 0.8 with a 1,000 bp window and a step size of 100 bp.
4. Remove all rows from the GWAS summary data that have duplicated rsIDs. After this, subset the quality-controlled UK Biobank genotype data and the GWAS data to a common set of variants, and flip the GWAS effect alleles to match ALT alleles in the UKB data. The fully processed genotype and GWAS data are saved to the subfolder `data_cleaned`.
5. Compute the top 20 genetic principal components from the quality-controlled genotype data. Following best practices, before computing PCs we further remove all variants in regions of long-range LD and then prune the remaining ones to a strict threshold of R^2 < 0.1 with a 1,000 bp window and a step size of 100 bp.
6. Create the file `data_cleaned/covariates.tsv` with one row per sample and 48 columns for sample ID, age, age^2, sex, age * sex, age^2 * sex, UKB assessment center (coded as 21 binary dummy variables), genotyping array (binary), and the first 20 genetic PCs. Finally, save the protein NPX levels to the file `data_cleaned/proteins.tsv` with one row per sample and one column per protein.

The `data_raw` folder can be deleted after the preprocessing script successfully finishes.

# Testing for association between co-expression and disease

This section describes how to perform a co-expression wide association study (COWAS), using data from the previous section as an example.

## Variant screening

For each protein, COWAS requires a list of genetic variants to be used as predictors in its expression imputation model. We considered three approaches for selecting variants:

1. **Variants selected by sure independence screening (SIS).** SIS is a variable selection method based on correlation learning. First, componentwise regression is performed to compute the marginal correlation between each standardized feature (genetic variant) and the response (protein expression). Then the features are ranked by their correlations, and the top $d$ (e.g. top 50) are selected as predictors.
2. **Variants that are pQTLs.** Instead of ranking features by their correlation with the response, they can be ranked by their marginal association *P*-value. That is, the top $d$ most significant pQTLs for the given protein can be used as predictors. Alternatively, one may wish to include all variants that pass a nominal significance threshold.
3. **Variants located near the gene coding for the given protein.** Since most of the genetic heritability of expression is explained by *cis*-QTLs, the two previous approaches can be restricted to variants that act locally on the gene coding for the given protein. For example, one may wish to include variants within a 1 Mb window of the gene boundaries. A protein annotation file listing, among other things, the start and end positions for all genes encoding the proteins included in UKB-PPP is available at <https://www.synapse.org/#!Synapse:syn52364558>. (Note that the positions in this file are on the GRCh38 build, while the UKB genotype data is on the GRCh37 build.)

We provide the script `map_pqtls.sh` for computing variant-protein associations for all proteins in the UKB-PPP dataset. The resulting summary statistics are saved to protein-specific, tab-separated files named `pqtls/<PROTEIN_CODE>_sumstats.tsv` with one line per variant and the following six columns: #CHROM, POS, ID, A1, BETA, P. These summary statistics can then be filtered to select predictors for the expression imputation models according to either the strength of their correlation (BETA) or the significance of their association (P).

Once you have determined which variants to use as predictors for each protein, save them to protein-specific files named `variants/<PROTEIN_CODE>_predictors.txt` with one rsID per line.

## Running COWAS

