# COWAS: The co-expression-wide association study

This repository provides a set of tools for conducting co-expression-wide association studies (COWAS). The goal of COWAS is to identify pairs of genes or proteins whose genetic component of co-expression is associated with complex traits. By considering the genetic regulation of both expression and co-expression, our method is able to boost power relative to standard TWAS and PWAS while also disentangling direct and interaction effects.

GitHub repo: https://github.com/mykmal/cowas  
Model weights: https://www.synapse.org/cowas  
Paper: https://doi.org/10.1101/2024.10.02.24314813

COWAS is run on one pair of genes or proteins at a time. First, three models are trained on an individual-level reference dataset. One model predicts the expression of the first gene/protein, another model predicts the expression of the second gene/protein, and the third model predicts the conditional co-expression of the two genes/proteins. Then the fitted model weights are used to impute expression and co-expression into GWAS data for any trait of interest. Finally, the outcome trait is jointly tested for association with the imputed expression and co-expression levels. We also provide trained imputation models for COWAS, enabling association testing to be performed using only GWAS summary statistics and a linkage disequilibrium (LD) reference panel.

## Installation

1. Download and unpack the COWAS repository from GitHub.
```bash
wget https://github.com/mykmal/cowas/archive/refs/heads/main.zip
unzip main.zip && rm main.zip
mv cowas-main cowas && cd cowas
```
2. Launch R and install the required packages optparse and data.table. If you wish to train your own imputation models then also install the package glmnet. If you wish to utilize parallel computation in glmnet then also install the package doMC. We used R 4.4.2, optparse 1.7.5, data.table 1.17.0, glmnet 4.1.8, and doMC 1.3.8.
```R
install.packages(c("optparse", "data.table", "glmnet", "doMC"))
```
3. Download PLINK 2.0 and place it in a directory on your PATH. We used PLINK v2.0.0-a.6.9LM AVX2 AMD (29 Jan 2025).
```bash
wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_amd_avx2_20250129.zip
unzip plink2_linux_amd_avx2_20250129.zip && rm plink2_linux_amd_avx2_20250129.zip
sudo mv plink2 /usr/local/bin/ && rm vcf_subset
```

## Usage guide

This section describes how to perform a co-expression-wide association study (COWAS) for any disease or phenotype of interest.

### Conducting association tests using COWAS

To perform a COWAS analysis, you will need to provide the following three ingredients:

1. Trained weights for expression and co-expression imputation models
2. GWAS summary statistics for your outcome trait of interest
3. Reference genotype data for computing a linkage disequilibrium (LD) panel

We provide weights for expression and co-expression imputation models trained on UK Biobank proteomic data at https://www.synapse.org/cowas. Several different sets of weights are available that were trained using different model choices or variant screening strategies, as described in our paper. Alternatively, you can train your own imputation models by following the instructions in the next section. The GWAS summary statistics file must contain the following columns: variant_id, effect_allele, other_allele, z_score, n_samples. The reference genotype data must be stored in a plain-text file with 0..2 coding. Importantly, the coding alleles in your genotype file must match the COWAS model effect alleles. If you are using our weights, this can be ensured by creating the reference genotype file using the `plink2 --export A` command with the `export-allele` option set to the ALT alleles in the `cowas_model_alleles.tsv` file we provide alongside our weights.

In addition to the above three ingredients, COWAS will ask for the name of each gene/protein in the current pair. You can also optionally set the output filename and the number of cores to use. Run `./cowas.R --help` to see the required syntax. The `cowas.R` script must be run separately for each pair of genes/proteins that you wish to analyze. To perform COWAS association testing for an entire list of pairs, use the batch script `run_cowas.sh` in the `utils` folder of this repository. Instructions are provided in comments at the top of the script.

Once COWAS finishes running, it will save association testing results to a tab-separated file with one line per pair. The output file columns have the following definitions:

| Column                 | Description |
| ---------------------- | ----------- |
| ID_A                   | Name of the first gene/protein in the pair |
| ID_B                   | Name of the second gene/protein in the pair |
| N_REFERENCE            | Sample size of the LD reference panel |
| N_GWAS                 | Sample size of the outcome trait GWAS |
| THETA_MARGINAL_A       | Marginal effect of the first gene/protein on the trait (identical to standard TWAS/PWAS) |
| VAR_THETA_MARGINAL_A   | Variance of the marginal effect size for the first gene/protein |
| PVAL_THETA_MARGINAL_A  | *P* value of the marginal effect of the first gene/protein |
| THETA_MARGINAL_B       | Marginal effect of the second gene/protein on the trait (identical to standard TWAS/PWAS) |
| VAR_THETA_MARGINAL_B   | Variance of the marginal effect size for the second gene/protein |
| PVAL_THETA_MARGINAL_B  | *P* value of the marginal effect of the second gene/protein |
| THETA_MARGINAL_CO      | Marginal effect of co-expression on the trait |
| VAR_THETA_MARGINAL_CO  | Variance of the marginal effect size of co-expression |
| PVAL_THETA_MARGINAL_CO | *P* value of the marginal effect of co-expression |
| THETA_JOINT_A          | Effect of the first gene/protein on the trait in a joint COWAS model |
| VAR_THETA_JOINT_A      | Variance of the first gene/protein effect size in a joint COWAS model |
| PVAL_THETA_JOINT_A     | *P* value of the first gene/protein in a joint COWAS model |
| THETA_JOINT_B          | Effect of the second gene/protein on the trait in a joint COWAS model |
| VAR_THETA_JOINT_B      | Variance of the second gene/protein effect size in a joint COWAS model |
| PVAL_THETA_JOINT_B     | *P* value of the second gene/protein in a joint COWAS model |
| THETA_JOINT_CO         | Effect of co-expression on the trait in a joint COWAS model |
| VAR_THETA_JOINT_CO     | Variance of the effect size of co-expression in a joint COWAS model |
| PVAL_THETA_JOINT_CO    | *P* value of co-expression in a joint COWAS model |
| FSTAT_JOINT            | *F* statistic for testing the overall significance of the joint COWAS model |
| PVAL_FSTAT_JOINT       | *P* value for overall significance of the joint COWAS model |

The columns named PVAL_FSTAT_JOINT and PVAL_THETA_JOINT_CO correspond to *P* values for the COWAS global test and the COWAS interaction test described in our paper, respectively. The *P* values in columns PVAL_THETA_JOINT_A and PVAL_THETA_JOINT_B can reveal whether each gene/protein is significant after accounting for the other gene/protein and their interaction, while the *P* values in columns PVAL_THETA_MARGINAL_A and PVAL_THETA_MARGINAL_B are equivalent to what one would get from a standard TWAS/PWAS analysis.

### Training your own imputation models

If you have access to a dataset with individual-level genotype data and gene or protein expression data, you can train your own models for imputing expression and co-expression. Model training is performed for one pair of genes/proteins at a time using the script `cowas_train.R`. To view the required inputs and the documentation for each command-line option, run `./cowas_train.R --help` from your command line.

The `cowas_train.R` script will save fitted model weights to an RDS file named `<ID_A>-<ID_B>.weights.rds` in the specified output folder. The RDS file stores a list of three named vectors, which contain genetic variant weights for the three models. In addition to saving model weights, `cowas_train.R` will also write model performance metrics to a tab-separated file named `performance_metrics.tsv` within the specified output directory. (Note that if this file already exists, a new line will be appended to its end.) The performance metrics file contains one line for each pair and the following nine columns:

| Column         | Description |
| -------------- | ----------- |
| ID_A           | Name or identifier of the first gene/protein |
| ID_B           | Name or identifier of the second gene/protein |
| SAMPLE_SIZE    | Sample size of the training data |
| NFEATURES_A    | Number of variants with nonzero weights in the <ID_A> model |
| CORRELATION_A  | Correlation between measured and predicted expression for <ID_A>, evaluated on a held-out 20% test set |
| NFEATURES_B    | Number of variants with nonzero weights in the <ID_B> model |
| CORRELATION_B  | Correlation between measured and predicted expression for <ID_B>, evaluated on a held-out 20% test set |
| NFEATURES_CO   | Number of variants with nonzero weights in the co-expression model |
| CORRELATION_CO | Correlation between estimated and predicted co-expression, evaluated on a held-out 20% test set |

We also provide the batch script `run_cowas_train.sh`, located in the `utils` folder of this repository, to automate the model training process for a user-specified list of gene or protein pairs. For each pair, this shell script will extract predictor variants from a PLINK-format file, convert their genotypes to the required format for COWAS, and then run `cowas_train.R`. Instructions are provided in comments at the top of the script.

### Variant screening tips

The `cowas_train.R` script will use all variants in the provided genotype matrices as model inputs. To reduce the computational time needed for model training, we suggest selecting an initial set of genetic variants for each gene or protein before running `cowas_train.R` or `run_cowas_train.sh`. In our paper, we compared three approaches for pre-screening variants:

1. **Variants selected by *P* value.** One approach is to select genetic variants that have the most significant association with expression levels. That is, the top $d$ (e.g., top 100) most significant QTLs for the given gene/protein can be used as an initial set of predictors. Alternatively, one may consider including all variants that pass a nominal significance threshold.
2. **Variants selected by effect size.** Instead of ranking variants by their QTL *P* values, they can be ranked by their marginal correlation with expression levels. First, componentwise regression is performed to compute the marginal correlation between each standardized feature (genetic variant) and the response (expression levels). Then the features are ranked by the absolute values of their correlations, and the top $d$ (e.g., top 100) are selected as predictors.
3. **Variants located in the *cis* region.** Since most of the genetic heritability of gene and protein expression is explained by *cis*-QTLs, the two previous approaches can be restricted to variants that act locally on the given gene or protein. In our paper, we considered variants within a 500 kb window of the gene boundaries for the gene coding for each protein. If you are using the UK Biobank plasma proteomics data, you can find the start and end positions for all genes encoding the assayed proteins at <https://www.synapse.org/#!Synapse:syn52364558>. (Note, however, that the positions given in that file are on the GRCh38 build, while positions in the UK Biobank genotype data are on the GRCh37 build.)

We provide the shell script `utils/map_pqtls.sh` for computing variant-protein associations for all proteins in the UK Biobank plasma proteomics dataset, and it can be easily modified for use with other datasets as well. The resulting summary statistics are saved to protein-specific, tab-separated files named `pqtls/<GENE_NAME>.sumstats.tsv` with one line per variant and the following six columns: #CHROM, POS, ID, A1, BETA, P. These summary statistics can then be filtered to select predictors for the expression imputation models according to either the strength of their correlation (BETA) or the significance of their association (P). For our paper, we pre-screened predictors using the scripts `utils/extract_predictors_general.R` and `utils/extract_predictors_specific.R`.

## Appendix: Data preparation and QC

This section describes how to obtain and prepare the data we used in our paper.

### UK Biobank downloads

To obtain individual-level data from the UK Biobank you will need to submit an application through the [Access Management System (AMS)](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access). After your application is approved, follow the UK Biobank documentation to download, unpack, and extract the following data fields:

* Data-Field 31 within Category 100094 (self-reported sex)
* Data-Field 54 within Category 100024 (UK Biobank assessment center)
* Data-Field 21003 within Category 100024 (age when attended assessment center)
* Data-Field 22000 within Category 100313 (genotype measurement batch)
* Data-Field 22006 within Category 100313 (indicator for individuals who self-identified as White British and have very similar genetic ancestry)
* Data-Field 22020 within Category 100313 (indicator for unrelated, high-quality samples used in PCA calculation)
* Data-Field 22828 within Category 100319 (WTCHG imputed genotypes)
* Data-Field 30900 within Category 1839 (UK Biobank plasma proteomics data)

Next, use PLINK to convert the genotype data to PLINK 2.0 binary format. Our QC scripts assume that the genotype data is stored in chromosome-specific files with filenames `ukb_chr<CHR>.pgen` + `ukb_chr<CHR>.psam` + `ukb_chr<CHR>.pvar`. Furthermore, we assume that the proteomic data is stored in a tab-separated file named `olink_data.tsv` beginning with a header line and containing the following four columns: eid, ins_index, protein_id, result. Finally, we assume that the remaining data fields listed above are stored in a single tab-separated file named `ukb_main_dataset.tsv` beginning with a header line and containing the following 13 columns: f.eid, f.31.0.0, f.54.0.0, f.54.1.0, f.54.2.0, f.54.3.0, f.21003.0.0, f.21003.1.0, f.21003.2.0, f.21003.3.0, f.22000.0.0, f.22006.0.0, f.22020.0.0. Copy all of these files to a new subfolder named `data_raw` within your main COWAS folder.

You will also need Data-Coding 143, which is a flat list containing gene names for each assayed protein. Download the file `coding143.tsv` from <https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=143> and place it in `data_raw`.

### GWAS summary statistics

For LDL cholesterol levels, we used summary statistics data from the Global Lipids Genetics Consortium (GLGC) GWAS conducted by [Graham et al. (2021)](https://www.nature.com/articles/s41586-021-04064-3). Note that this study provides multi-ancestry as well as ancestry-specific results, but we only considered the European results in order to match the genetic ancestry of the UK Biobank. The summary-level associations can be downloaded from <https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/>. Place the unpacked file `LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results` in `data_raw`.

For our primary analysis of Alzheimer's disease, we used summary statistics data from the European Alzheimer & Dementia Biobank (EADB) consortium GWAS conducted by [Bellenguez et al. (2022)](https://www.nature.com/articles/s41588-022-01024-z). This was the largest genetic association study of Alzheimer's disease available at the time of writing. Summary-level associations for the EADB GWAS can be downloaded from <https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/>. Place the unpacked file `GCST90027158_buildGRCh38.tsv` in `data_raw`.

Although the EADB GWAS of Alzheimer's disease is the largest available, some have questioned its quality because the EADB study relied in part on proxy cases, where participants' disease status was imputed from their family history of Alzheimer's disease instead of being clinically diagnosed. Thus, we also considered the largest GWAS of Alzheimer's disease that does not contain proxy cases: the International Genomics of Alzheimer's Project (IGAP) consortium GWAS, which was coducted by [Kunkle et al. (2019)](https://www.nature.com/articles/s41588-019-0358-2). Summary-level associations for the IGAP GWAS can be downloaded from <https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007511>. Place the file `Kunkle_etal_Stage1_results.txt` in `data_raw`.

For Parkinson's disease, we used summary statistics data from the International Parkinson Disease Genomics Consortium (IPDGC) GWAS conducted by [Nalls et al. (2019)](https://www.thelancet.com/journals/laneur/article/PIIS1474-4422(19)30320-5). The summary-level associations can be downloaded from <https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/harmonised/>. Place the unpacked file `GCST009325.h.tsv` in `data_raw`.

### Data processing

Run the shell script `qc/preprocess_ukb.sh` from within the main COWAS folder to process the downloaded data. This script will perform the following data wrangling and quality control steps:

1. Create the file `pairs/all_protein_pairs.tsv` listing all possible pairs of proteins that have data available, with one pair per row.
2. Filter the main dataset to obtain a set of high-quality, unrelated, White British samples with per-sample genotyping rate > 99%. Filter the proteomic data to obtain the set of samples assessed at the initial visit. Then subset the genotype data, proteomic data, and covariate data to a common set of samples.
3. Remove variants from the UK Biobank genotype data that have a missingness rate $> 10\%$, have an MAC $< 100$, have an MAF $< 1\%$, fail a Hardy-Weinberg equilibrium test (*P* $< 10^{-15}$), lack an rsID, or are palindromic. Next, prune the remaining variants to $R^2 < 0.8$ with a 1,000 bp window and a step size of 100 bp. The quality-controlled genotype data is saved to the files `genotypes.pgen` + `genotypes.psam` + `genotypes.pvar` within the subfolder `data_cleaned`.
4. Compute the top 20 genetic principal components from the quality-controlled genotype data. Following best practices, before computing PCs we further remove all variants in regions of long-range LD and then prune the remaining ones to a strict threshold of $R^2 < 0.1$ with a 1,000 bp window and a step size of 100 bp.
5. Create the file `data_cleaned/covariates.tsv` with one row per sample and 48 columns for sample ID, age, age^2, sex, age \* sex, age^2 \* sex, UK Biobank assessment center (coded as 21 binary dummy variables), genotyping array (binary), and the first 20 genetic PCs. Save the protein NPX levels to the file `data_cleaned/proteins.tsv` with one row per sample and one column per protein.
6. For each GWAS, remove any rows that have duplicated rsIDs. After this, subset the quality-controlled genotype data and the GWAS data to a common set of variants, and flip the GWAS effect alleles and effect sizes to match ALT alleles in the genotype data. The fully processed GWAS summary statistics, as well as copies of the genotype data subset to variants present in each GWAS, are saved to the subfolder `data_cleaned`.

The `data_raw` subfolder can be deleted after the preprocessing script successfully finishes.

