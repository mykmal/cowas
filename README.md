# Overview of COWAS

This repository provides a set of tools for performing co-expression wide association studies (COWAS).

The goal of COWAS is to test for association between disease and the genetically regulated componenet of gene co-expression. First, models are trained to predict gene expression and co-expression from genotype data. Next, the fitted model weights are used to impute gene expression and co-expression into disease GWAS data. Finally, disease is jointly tested for association with the expression level of each gene as well as with their co-expression. COWAS can be conceptually thought of as fitting the following regression model, where the interaction term represents gene co-expression:
```
disease ~ gene_1_expression + gene_2_expression + gene_1_expression * gene_2_expression
```
If the coefficient for the interaction term is significant, we conclude that genetically regulated co-expression between the given pair of genes influences disease risk. The direction and magnitude of this effect are given by the sign and magnitude of the coefficent estimate, respectively.

# Software documentation

## Installation

1. Download and unpack the COWAS repository from GitHub.
```bash
wget https://github.com/MykMal/cowas/archive/refs/heads/main.zip
unzip main.zip && rm main.zip
mv cowas-main cowas && cd cowas
```
2. Launch R and install the packages optparse, data.table, and glmnet.
```R
install.packages(c("optparse", "data.table", "glmnet"))
```
3. If you wish to train your own co-expression imputation models, download PLINK 1.90 to the unpacked COWAS folder. We used PLINK v1.90b7.1 64-bit (18 Oct 2023).
```bash
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231018.zip
unzip plink_linux_x86_64_20231018.zip plink
rm plink_linux_x86_64_20231018.zip
```

## Running COWAS with pre-computed models

## Training co-expression imputation models

# Appendix: GTEx data pre-processing

This appendix describes how to obtain and prepare the data used in our paper.

From within the unpacked COWAS folder, run the following terminal commands to download GTEx v8 expression data, covariate information, and gene annotations.
```bash
mkdir raw && cd raw
wget https://storage.cloud.google.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_EUR.tar
tar -xf GTEx_Analysis_v8_eQTL_EUR.tar
wget https://storage.cloud.google.com/adult-gtex/references/v8/gencode.v26.GRCh38.genes.gtf
```

To obtain individual-level genotype data from the GTEx Project you will need to submit an application through [dbGaP](https://www.ncbi.nlm.nih.gov/gap/). After obtaining access to the data under accession number phs000424.v8.p2, follow the dbGaP documentation to download the following files and extract their contents to `cowas/raw`:

* `phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar'
* `phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz`

Finally, run the shell script `util/preprocess_gtex.sh` from within the main COWAS folder. This will perform extensive quality control steps on the genotype data, reformat the expression matrices and expression covariates, create an annotation file with all GTEx gene transcripts, and create files listing all gene pairs available for each tissue. The `cowas/raw` folder can be deleted after the script successfully finishes.


