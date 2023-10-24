#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --time=96:00:00
#SBATCH --partition=msismall,msilarge,msibigmem
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out

module load R/4.3.0-openblas

# Specify a file that lists pairs of genes for which to train COWAS weights.
# This should be a text file with one tab-separated pair of gene names per line.
GENES=pairs/Whole_Blood_chr1_part000

# Suffix that will be appended to the name of the .pos output file.
# This is useful when a single chromosome is split across several gene pair lists.
SUFFIX=000

# Chromosome number corresponding to the gene pairs in $GENES
CHR=1

# Name of the tissue for which to train models.
# This should match the prefix of the relevant expression/covariate file names.
TISSUE=Whole_Blood

# Directory for storing COWAS weights and performance metrics.
OUT=/scratch.global/malak039_output


printf "Runtime parameters:\n"
printf "GENES = ${GENES}\n"
printf "SUFFIX = ${SUFFIX}\n"
printf "CHR = ${CHR}\n"
printf "TISSUE = ${TISSUE}\n"
printf "OUT = ${OUT}\n\n"

if ( [ ! -f ${OUT}/GTExv8.EUR.${TISSUE}.CHR${CHR}_${SUFFIX}.pos ] ); then
printf "PANEL\tWGT\tID\tCHR\tP0\tP1\tN\n" > ${OUT}/GTExv8.EUR.${TISSUE}.CHR${CHR}_${SUFFIX}.pos
fi

# Reads names of genes with available expression measurements
read -r MEASURED_GENES < expression/${TISSUE}.expression.txt

# Loops through all gene pairs in the $GENES file
while read -r GENE_A GENE_B ETC; do

if ( [[ ${MEASURED_GENES} != *"${GENE_A}"* ]] || [[ ${MEASURED_GENES} != *"${GENE_B}"* ]] ); then
printf "WARNING: expression data not found for ${GENE_A} or ${GENE_B}. Skipping this pair.\n"
continue
fi

# Extracts gene info from annotation file
read CHR_A START_A END_A <<<$(awk -v OFS='\t' -v ga=${GENE_A} '$2 == ga {print $3,$4,$5}' annotations.txt)
read CHR_B START_B END_B <<<$(awk -v OFS='\t' -v gb=${GENE_B} '$2 == gb {print $3,$4,$5}' annotations.txt)

# Removes chr prefix (if it exists)
CHR_A=$(echo ${CHR_A} | tr -d "chr")
CHR_B=$(echo ${CHR_B} | tr -d "chr")

if ( [[ ${CHR_A} != ${CHR} ]] || [[ ${CHR_B} != ${CHR} ]] ); then
printf "WARNING: ${GENE_A} or ${GENE_B} is not located on chromosome ${CHR}. Skipping this pair.\n"
continue
fi

# Creates folder for storing temporary files
if [ ! -d "${OUT}/temp" ]; then
mkdir ${OUT}/temp
fi
COWAS_TEMP=${OUT}/temp/${GENE_A}_${GENE_B}
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

# Extracts the cis-region, subsets to individuals with expression values,
# then creates a bim file and exports genotypes to a text file with 0/1/2 coding
./plink --bfile genotypes/dosages_chr${CHR} \
          --silent \
          --chr ${CHR} \
          --from-bp ${START_WINDOW_A} \
          --to-bp ${END_WINDOW_A} \
          --keep expression/${TISSUE}.individuals.txt \
          --recode 01 A tabx \
          --make-just-bim \
          --out ${COWAS_TEMP}/${GENE_A}
./plink --bfile genotypes/dosages_chr${CHR} \
          --silent \
          --chr ${CHR} \
          --from-bp ${START_WINDOW_B} \
          --to-bp ${END_WINDOW_B} \
          --keep expression/${TISSUE}.individuals.txt \
          --recode 01 A tabx \
          --make-just-bim \
          --out ${COWAS_TEMP}/${GENE_B}

if ( [ ! -f ${COWAS_TEMP}/${GENE_A}.raw ] || [ ! -f ${COWAS_TEMP}/${GENE_B}.raw ] ); then
printf "WARNING: Unable to extract genotype data for ${GENE_A} or ${GENE_B}. Skipping this pair.\n"
rm -rf ${COWAS_TEMP}
continue
fi

# Removes unnecessary columns from the raw PLINK genotype files
cut -f 2,7- ${COWAS_TEMP}/${GENE_A}.raw > ${COWAS_TEMP}/${GENE_A}.gmatrix
cut -f 2,7- ${COWAS_TEMP}/${GENE_B}.raw > ${COWAS_TEMP}/${GENE_B}.gmatrix

# Creates a single bim file with variants for both genes.
# This will duplicate SNPs in overlapping cis-regions, but cowas.R takes care of that later.
cat ${COWAS_TEMP}/${GENE_A}.bim ${COWAS_TEMP}/${GENE_B}.bim > ${COWAS_TEMP}/both.snps

./cowas_train.R --gene_a ${GENE_A} \
          --gene_b ${GENE_B} \
          --genotypes_a ${COWAS_TEMP}/${GENE_A}.gmatrix \
          --genotypes_b ${COWAS_TEMP}/${GENE_B}.gmatrix \
          --expression expression/${TISSUE}.expression.txt \
          --covariates covariates/${TISSUE}.covariates.txt \
          --snps ${COWAS_TEMP}/both.snps \
          --out ${OUT}/GTExv8.EUR.${TISSUE}.CHR${CHR}_${SUFFIX}

if ( [ -f ${OUT}/GTExv8.EUR.${TISSUE}.CHR${CHR}_${SUFFIX}/${GENE_A}_${GENE_B}_wgt.RData ] ); then
N=$( wc -l expression/${TISSUE}.individuals.txt | awk '{print $1}' )
printf "GTExv8.EUR.${TISSUE}.CHR${CHR}\tGTExv8.EUR.${TISSUE}.CHR${CHR}/${GENE_A}_${GENE_B}_wgt.RData\t${GENE_A}_${GENE_B}\t${CHR}\t${START_A}\t${START_B}\t${N}\n" \
          >> ${OUT}/GTExv8.EUR.${TISSUE}.CHR${CHR}_${SUFFIX}.pos
fi

rm -rf ${COWAS_TEMP}

done < ${GENES}

printf "Done!\n"

