#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 10
#SBATCH --job-name=g2gtools
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/g2gtools.txt
#SBATCH -e logs/g2gtools.txt
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"

ml load singularity 
# Define the path to the g2gtools Singularity image
SIF_PATH=/users/sparthib/retina_lrs/code/08_ASE/g2g_tools/g2gtools.sif

# Define input directories and files
VCF_INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf
SNP=${VCF_INPUT_DIR}/filtered_SNP.vcf
INDEL=${VCF_INPUT_DIR}/filtered_INDEL.vcf

REF_DIR=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly
FASTA=${REF_DIR}/release_46_primary_genome.fasta
GTF=${REF_DIR}/release_46_primary_assembly.gtf

OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/g2gtools

# Define strain names and output directories
unfiltered_vcf=${VCF_INPUT_DIR}/genotyped.vcf.gz
STRAIN1_NAME=SRR1091088
STRAIN2_NAME=SRR1091091
STRAIN1_DIR=$OUTPUT_DIR/SRR1091088
STRAIN2_DIR=$OUTPUT_DIR/SRR1091091

# Create output directories
# mkdir -p $STRAIN1_DIR
# mkdir -p $STRAIN2_DIR

# VCF to VCI conversion
echo "convert vcf to vci"
singularity exec $SIF_PATH g2gtools vcf2vci --help

# Run the actual vcf2vci command with proper directory bindings
# singularity exec --bind ${REF_DIR}:${REF_DIR},${VCF_INPUT_DIR}:${VCF_INPUT_DIR},${STRAIN1_DIR}:${STRAIN1_DIR} $SIF_PATH \
# g2gtools vcf2vci -i $unfiltered_vcf -f $FASTA -s ${STRAIN1_NAME} -o ${STRAIN1_DIR}/output.vci --diploid --pass

grep -v '*' ${STRAIN1_DIR}/output.vci  > ${STRAIN1_DIR}/output_filtered.vci #remove lines that have * allele
gzip ${STRAIN1_DIR}/output_filtered.vci #compress the file

output_filtered.vci.gz
# Patching the genome with SNPs from the VCI
echo "patching genome"
singularity exec --bind ${REF_DIR}:${REF_DIR},$STRAIN1_DIR:$STRAIN1_DIR $SIF_PATH \
g2gtools patch -i $FASTA -c ${STRAIN1_DIR}/output_filtered.vci.gz -o ${STRAIN1_DIR}/${STRAIN1_NAME}_patched.fa

# Transforming the patched genome to diploid
echo "transforming genome to diploid"
singularity exec --bind ${STRAIN1_DIR}:${STRAIN1_DIR} $SIF_PATH \
g2gtools transform -i ${STRAIN1_DIR}/${STRAIN1_NAME}_patched.fa -c ${STRAIN1_DIR}/output_filtered.vci.gz -o ${STRAIN1_DIR}/${STRAIN1_NAME}_diploid_genome.fa

# Converting the VCI to GTF format
echo "convert vci to gtf"
singularity exec --bind ${STRAIN1_DIR}:${STRAIN1_DIR},${REF_DIR}:${REF_DIR} $SIF_PATH \
g2gtools convert -c ${STRAIN1_DIR}/output_filtered.vci.gz -i $GTF -f GTF -o ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf

# Creating a GTF database
echo "create gtf database"
singularity exec --bind ${STRAIN1_DIR}:${STRAIN1_DIR} $SIF_PATH \
g2gtools gtf2db -i ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf -o ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf.db


# Extracting transcripts from the diploid genome
echo "extract transcripts"
singularity exec --bind ${STRAIN1_DIR}:${STRAIN1_DIR} $SIF_PATH \
g2gtools extract --transcripts -i ${STRAIN1_DIR}/${STRAIN1_NAME}_diploid_genome.fa -db ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf.db > ${STRAIN1_DIR}/${STRAIN1_NAME}.transcripts.fa


# use prepare-emase function  from emase package (not in g2gtools)
# to create diploid transcriptome from NOD and PWk genome and GTF files
# GENOME1=${STRAIN1}/${STRAIN1}.fa
# GENOME2=${STRAIN2}/${STRAIN2}.fa
# GTF1=${STRAIN1}/${STRAIN1}.gtf
# GTF2=${STRAIN2}/${STRAIN2}.gtf
# SUFFIX1=N
# SUFFIX2=P
# EMASE_DIR=NxP
# prepare-emase -G ${GENOME1},${GENOME2} -g ${GTF1},${GTF2} -s ${SUFFIX1},${SUFFIX2} -o ${EMASE_DIR} -m

echo "**** Job ends ****"
date +"%Y-%m-%d %T"