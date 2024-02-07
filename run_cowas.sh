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

# A file that lists pairs of genes or proteins for which to train COWAS weights
# This should be a text file with one tab-separated pair of gene/protein names per line
PAIRS=phenotypes/protein_pairs.tsv

# File name for COWAS results
OUT_FILE=ukb_proteins_all.tsv

# Directory for storing the COWAS output file
OUT_DIR=output


printf "Runtime parameters:\n"
printf "PAIRS = ${PAIRS}\n"
printf "OUT_FILE = ${OUT_FILE}\n"
printf "OUT_DIR = ${OUT_DIR}\n\n"

if ( [ ! -d ${OUT_DIR} ] ); then
mkdir ${OUT_DIR}
fi

if ( [ ! -f ${OUT_DIR}/${OUT_FILE} ] ); then
printf "ID_A\tID_B\tN_REFERENCE\tN_GWAS\t\
NFEATURES_A\tR2PRED_A\t\
NFEATURES_B\tR2PRED_B\t\
NFEATURES_CO\tR2PRED_CO\t\
THETA_DIRECT_A\tVAR_THETA_DIRECT_A\tPVAL_THETA_DIRECT_A\t\
THETA_DIRECT_B\tVAR_THETA_DIRECT_B\tPVAL_THETA_DIRECT_B\t\
THETA_DIRECT_CO\tVAR_THETA_DIRECT_CO\tPVAL_THETA_DIRECT_CO\t\
THETA_FULL_A\tVAR_THETA_FULL_A\tPVAL_THETA_FULL_A\t\
THETA_FULL_B\tVAR_THETA_FULL_B\tPVAL_THETA_FULL_B\t\
THETA_FULL_CO\tVAR_THETA_FULL_CO\tPVAL_THETA_FULL_CO\n" > ${OUT_DIR}/${OUT_FILE}
fi

# Reads names of proteins with available expression measurements
read -r MEASURED_PROTEINS < phenotypes/proteins.tsv

# Loops through all protein pairs in the $PAIRS file
while read -r PROTEIN_A PROTEIN_B ETC; do

if ( [[ ${PROTEIN_A} == ${PROTEIN_B} ]] ); then
continue
fi

if ( [[ ${MEASURED_PROTEINS} != *"${PROTEIN_A}"* ]] || [[ ${MEASURED_PROTEINS} != *"${PROTEIN_B}"* ]] ); then
printf "WARNING: expression data not found for ${PROTEIN_A} or ${PROTEIN_B}. Skipping this pair.\n"
continue
fi

# Extracts protein info from the annotation file
CODE_A=$(echo ${PROTEIN_A} | tr -d "protein_")
CODE_B=$(echo ${PROTEIN_B} | tr -d "protein_")
read CHR_A START_A END_A <<<$(awk -v FS='\t' -v OFS='\t' -v pa=${CODE_A} '(NR > 1) && ($1 == pa) {print $5,$6,$7}' util/olink_annotations.tsv)
read CHR_B START_B END_B <<<$(awk -v FS='\t' -v OFS='\t' -v pb=${CODE_B} '(NR > 1) && ($1 == pb) {print $5,$6,$7}' util/olink_annotations.tsv)

# Removes chr prefix (if it exists)
CHR_A=$(echo ${CHR_A} | tr -d "chr")
CHR_B=$(echo ${CHR_B} | tr -d "chr")

# Creates folder for storing temporary files
if [ ! -d "${OUT_DIR}/temp" ]; then
mkdir ${OUT_DIR}/temp
fi
COWAS_TEMP=${OUT_DIR}/temp/${CODE_A}_${CODE_B}
mkdir ${COWAS_TEMP}

# Expands 1 Mb before the transcription start site
START_WINDOW_A=$(( ${START_A} - 1000000 ))
if (( ${START_WINDOW_A} < 0 )); then
START_WINDOW_A=0
fi
START_WINDOW_B=$(( ${START_B} - 1000000 ))
if (( ${START_WINDOW_B} < 0 )); then
START_WINDOW_B=0
fi

# Expands 1 Mb after the transcription end site
END_WINDOW_A=$(( ${END_A} + 1000000 ))
END_WINDOW_B=$(( ${END_B} + 1000000 ))

awk -v FS="\t" -v OFS="\t" '!/^#/ {print $3,$5}' genotypes/ukb_filtered_${CHR_A}.pvar > ${COWAS_TEMP}/A_alleles.tsv
awk -v FS="\t" -v OFS="\t" '!/^#/ {print $3,$5}' genotypes/ukb_filtered_${CHR_B}.pvar > ${COWAS_TEMP}/B_alleles.tsv

# Extracts the cis-region, subsets to individuals with expression values,
# then creates a bim file and exports genotypes to a text file with 0/1/2 coding
plink2 --pfile genotypes/ukb_filtered_${CHR_A} \
       --silent \
       --chr ${CHR_A} \
       --from-bp ${START_WINDOW_A} \
       --to-bp ${END_WINDOW_A} \
       --export A \
       --export-allele ${COWAS_TEMP}/A_alleles.tsv \
       --make-just-pvar cols=maybecm \
       --out ${COWAS_TEMP}/${PROTEIN_A}
plink2 --pfile genotypes/ukb_filtered_${CHR_B} \
       --silent \
       --chr ${CHR_B} \
       --from-bp ${START_WINDOW_B} \
       --to-bp ${END_WINDOW_B} \
       --export A \
       --export-allele ${COWAS_TEMP}/B_alleles.tsv \
       --make-just-pvar cols=maybecm \
       --out ${COWAS_TEMP}/${PROTEIN_B}

if ( [ ! -f ${COWAS_TEMP}/${PROTEIN_A}.raw ] || [ ! -f ${COWAS_TEMP}/${PROTEIN_B}.raw ] ); then
printf "WARNING: Unable to extract genotype data for ${PROTEIN_A} or ${PROTEIN_B}. Skipping this pair.\n"
rm -rf ${COWAS_TEMP}
continue
fi

# Removes unnecessary columns from the raw PLINK genotype files
cut -f 2,7- ${COWAS_TEMP}/${PROTEIN_A}.raw > ${COWAS_TEMP}/${PROTEIN_A}.gmatrix
cut -f 2,7- ${COWAS_TEMP}/${PROTEIN_B}.raw > ${COWAS_TEMP}/${PROTEIN_B}.gmatrix

# Creates a single pvar file with variants for both proteins
# This will duplicate SNPs in overlapping cis-regions, but cowas.R takes care of that later
cat ${COWAS_TEMP}/${PROTEIN_A}.pvar ${COWAS_TEMP}/${PROTEIN_B}.pvar > ${COWAS_TEMP}/both.snps

# And this is where the magic happens!
./cowas.R --protein_a ${PROTEIN_A} \
          --protein_b ${PROTEIN_B} \
          --genotypes_a ${COWAS_TEMP}/${PROTEIN_A}.gmatrix \
          --genotypes_b ${COWAS_TEMP}/${PROTEIN_B}.gmatrix \
          --snps ${COWAS_TEMP}/both.snps \
          --expression phenotypes/proteins.tsv \
          --covariates phenotypes/covariates.tsv \
          --gwas gwas/Bellenguez_2022_AD_gwas.tsv \
          --out ${OUT_DIR}/${OUT_FILE}

rm -rf ${COWAS_TEMP}

done < ${PAIRS}

printf "Done!\n"

