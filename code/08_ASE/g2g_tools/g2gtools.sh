#!/bin/bash

source activate g2gtools


# adapted from https://github.com/narayananr/diploid_txome/blob/master/create_diploid_transcriptome.sh

# get SNP and INDEL VCF files 
VCF_INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf
SNP=${VCF_INPUT_DIR}/filtered_SNP.vcf
INDEL=${VCF_INPUT_DIR}/filtered_INDEL.vcf
GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf
REF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_all_transcripts.fa
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
g2gtools vcf2vci -f ${REF} -i $unfiltered_vcf -s ${STRAIN1_NAME} -o ${STRAIN1_DIR}/output.vci --diploid
g2gtools patch -i $REF -c ${STRAIN1_DIR}/output.vci -o ${STRAIN1_DIR}/${STRAIN1_NAME}_patched.fa
g2gtools transform -i ${STRAIN1_DIR}/${STRAIN1_NAME}_patched.fa -c ${STRAIN1_DIR}/output.vci -o ${STRAIN1_DIR}/${STRAIN1_NAME}_diploid_genome.fa


g2gtools convert -c ${STRAIN1_DIR}/output.vci -i ${GTF}  -o ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf
g2gtools gtf2db -i ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf -o ${STRAIN1_DIR}/${STRAIN1_NAME}.gtf.db

# extract transcripts from NOD genome
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