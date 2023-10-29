# Overview of COWAS

This repository provides a set of tools for performing co-expression wide association studies (COWAS).

The goal of COWAS is to test for association between disease and the genetically regulated component of gene co-expression. First, models are trained to predict gene expression and co-expression from genotype data. Then the fitted model weights are used to impute gene expression and co-expression into disease GWAS data. Finally, disease is jointly tested for association with the expression level of each gene as well as with their co-expression. By considering the effect of correlated gene expression on disease in addition to direct effects, COWAS can identify novel disease-relevant genes and aid in the interpretation of GWAS findings.

# Software documentation

## Installation

1. Download and unpack the COWAS repository from GitHub.
```bash
wget https://github.com/MykMal/cowas/archive/refs/heads/main.zip
unzip cowas-main.zip && rm cowas-main.zip
mv cowas-main cowas && cd cowas
```
2. Launch R and install the packages optparse, data.table, and glmnet. We used R 4.3.0, optparse 1.7.3, data.table 1.14.8, and glmnet 4.1.8.
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

First, create a subfolder named `raw` within your unpacked COWAS folder. Then download the following files from the [GTEx Portal](https://www.gtexportal.org/home/downloads/adult-gtex) to `cowas/raw`:

* `GTEx_Analysis_v8_eQTL_EUR.tar` (from the "Single-Tissue cis-QTL Data" section of the QTL tab)
* `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz` (from the Reference tab)
* `gencode.v26.GRCh38.genes.gtf` (from the Reference tab)

To obtain individual-level genotype data from the GTEx Project you will need to submit an application through [dbGaP](https://www.ncbi.nlm.nih.gov/gap/). After obtaining access to the data under accession number phs000424.v8.p2, follow the dbGaP documentation to download and extract the following file to `cowas/raw`:

* `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz` (from the archive `phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar')

Finally, run the shell script `util/preprocess_gtex.sh` from within the main COWAS folder. This will perform extensive quality control steps on the genotype data, reformat the expression matrices and expression covariates, create an annotation file with all GTEx gene transcripts, and create files listing all gene pairs available for each tissue. The `cowas/raw` folder can be deleted after the script successfully finishes.

