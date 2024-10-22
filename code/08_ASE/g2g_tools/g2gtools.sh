#!/bin/bash

source activate g2gtools


#https://github.com/narayananr/diploid_txome/blob/master/create_diploid_transcriptome.sh

# get SNP and INDEL VCF files 
VCF_INPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf
SNP=${VCF_INPUT_DIR}/filtered_SNP.vcf
INDEL=${VCF_INPUT_DIR}/filtered_INDEL.vcf

# reference genome in fasta format
# ftp://ftp.ensembl.org/pub/release-75/gtf/mus_musculus/
REF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_all_transcripts.fa
# strain name (usually a column name in the vcf file), e.g., CAST_EiJ
STRAIN1=SRR1091088
STRAIN2=SRR1091091
mkdir -p ${SRR1091088}
mkdir -p ${SRR1091091}

###########
# g2gtools to create NOD specific genome and transcriptome
###########
# ftp://ftp.ensembl.org/pub/release-75/gtf/mus_musculus/

GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf
# Create a chain file for mapping bases between two genomes. In this case, between reference and NOD
g2gtools vcf2chain -f ${REF} -i ${NOD_PWK_Indel}.gz -s ${STRAIN1} -o ${STRAIN1}/REF-to-${STRAIN1}.chain
# patch SNPs on to reference genome
g2gtools patch -i ${REF} -s ${STRAIN1} -v ${NOD_PWK_SNP}.gz -o ${STRAIN1}/${STRAIN1}.patched.fa
g2gtools transform -i ${STRAIN1}/${STRAIN1}.patched.fa -c ${STRAIN1}/REF-to-${STRAIN1}.chain -o ${STRAIN1}/${STRAIN1}.fa
g2gtools convert -c ${STRAIN1}/REF-to-${STRAIN1}.chain -i ${GTF} -f gtf -o ${STRAIN1}/${STRAIN1}.gtf
#g2gtools gtf2db -i ${STRAIN1}/${STRAIN1}.gtf -o ${STRAIN1}/${STRAIN1}.gtf.db
# extract transcripts from NOD genome
#g2gtools extract --transcripts -i ${STRAIN1}/${STRAIN1}.fa -db ${STRAIN1}/${STRAIN1}.gtf.db > ${STRAIN1}/${STRAIN1}.transcripts.fa

###########
# g2gtools to create PWK specific genome and transcriptome
###########
# Create a chain file for mapping bases between two genomes. In this case, between reference and PWK
g2gtools vcf2chain -f ${REF} -i ${NOD_PWK_Indel}.gz -s ${STRAIN2} -o ${STRAIN2}/REF-to-${STRAIN2}.chain
# patch SNPs on to reference genome
g2gtools patch -i ${REF} -s ${STRAIN2} -v ${NOD_PWK_SNP}.gz -o ${STRAIN2}/${STRAIN2}.patched.fa
g2gtools transform -i ${STRAIN2}/${STRAIN2}.patched.fa -c ${STRAIN2}/REF-to-${STRAIN2}.chain -o ${STRAIN2}/${STRAIN2}.fa
# PWK specific GTF file
g2gtools convert -c ${STRAIN2}/REF-to-${STRAIN2}.chain -i ${GTF} -f gtf -o ${STRAIN2}/${STRAIN2}.gtf
#g2gtools gtf2db -i ${STRAIN2}/${STRAIN2}.gtf -o ${STRAIN2}/${STRAIN2}.gtf.db
# extract transcripts from PWK genome
#g2gtools extract --transcripts -i ${STRAIN2}/${STRAIN2}.fa -db ${STRAIN2}/${STRAIN2}.gtf.db > ${STRAIN2}/${STRAIN2}.transcripts.fa

# use prepare-emase function  from emase package (not in g2gtools)
# to create diploid transcriptome from NOD and PWk genome and GTF files
GENOME1=${STRAIN1}/${STRAIN1}.fa
GENOME2=${STRAIN2}/${STRAIN2}.fa
GTF1=${STRAIN1}/${STRAIN1}.gtf
GTF2=${STRAIN2}/${STRAIN2}.gtf
SUFFIX1=N
SUFFIX2=P
EMASE_DIR=NxP
prepare-emase -G ${GENOME1},${GENOME2} -g ${GTF1},${GTF2} -s ${SUFFIX1},${SUFFIX2} -o ${EMASE_DIR} -m