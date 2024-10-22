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
SIF_PATH=/users/sparthib/retina_lrs/code/08_ASE/g2g_tools/g2gtools.sif 

# adapted from https://github.com/narayananr/diploid_txome/blob/master/create_diploid_transcriptome.sh

# get SNP and INDEL VCF files 
VCF_INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf
SNP=${VCF_INPUT_DIR}/filtered_SNP.vcf
INDEL=${VCF_INPUT_DIR}/filtered_INDEL.vcf
GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf
REF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fasta
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/g2gtools

# strain name (usually a column name in the vcf file), e.g., CAST_EiJ
unfiltered_vcf=${VCF_INPUT_DIR}/genotyped.vcf
STRAIN1_NAME=SRR1091088
STRAIN2_NAME=SRR1091091
STRAIN1_DIR=$OUTPUT_DIR/SRR1091088
STRAIN2_DIR=$OUTPUT_DIR/SRR1091091
mkdir -p $STRAIN1_DIR
mkdir -p $STRAIN2_DIR

# Create a chain file for mapping bases between two genomes. In this case, between reference and NOD
# g2gtools vcf2chain -f ${REF} -i $INDEL -s ${STRAIN1_NAME} -o ${STRAIN1_DIR}/REF-to-${STRAIN1_NAME}.chain
# 
# # patch SNPs on to reference genome
# g2gtools patch -i ${REF} -s ${STRAIN1_NAME} -v $SNP -o ${STRAIN1_DIR}/${STRAIN1_NAME}.patched.fa
# g2gtools transform -i ${STRAIN1_DIR}/${STRAIN1_NAME}.patched.fa -c ${STRAIN1_DIR}/REF-to-${STRAIN1_NAME}.chain -o ${STRAIN1_DIR}/${STRAIN1_NAME}.fa


## try using vci instead 

echo "convert vcf to vci"
singularity exec $SIF_PATH g2gtools vcf2vci -f ${REF} -i $unfiltered_vcf -s ${STRAIN1_NAME} -o ${STRAIN1_DIR}/output.vci

echo "patching genome"
singularity exec $SIF_PATH g2gtools patch -i $REF -c ${STRAIN1_DIR}/output.vci -o ${STRAIN1_DIR}/${STRAIN1_NAME}_patched.fa

echo "transforming genome to diploid"
singularity exec $SIF_PATH g2gtools transform -i ${STRAIN1_DIR}/${STRAIN1_NAME}_patched.fa -c ${STRAIN1_DIR}/output.vci -o ${STRAIN1_DIR}/${STRAIN1_NAME}_diploid_genome.fa

echo "convert vci to gtf"
singularity exec $SIF_PATH g2gtools convert -c ${STRAIN1_DIR}/output.vci -i ${GTF} -o ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf

echo "create gtf database"
singularity exec $SIF_PATH g2gtools gtf2db -i ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf -o ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf.db

# extract transcripts from NOD genome
echo "extract transcripts"
singularity exec $SIF_PATH g2gtools extract --transcripts -i ${STRAIN1_DIR}/${STRAIN1_NAME}_diploid_genome.fa -db ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf.db > ${STRAIN1_DIR}/${STRAIN1_NAME}.transcripts.fa

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