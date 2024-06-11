#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 15
#SBATCH --job-name=sqanti
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/sqanti.txt
#SBATCH -e logs/sqanti.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"


#reference genome 
REFERENCE_GENOME_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
#reference transcriptome 
REFERENCE_TRANSCRIPTOME=/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts_short_header.fa

#sample transcriptome
ROs_GTF=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation/extended_annotation.gtf
FT_vs_RGC_GTF=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/extended_annotation.gtf
# GTF (default): by default, SQANTI3 expects the transcriptome to be provided as a GTF file, 
# and we recommend to stick to this format if your transcriptome construction pipeline allows it
RO_OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation/sqanti3_qc
FT_vs_RGC_OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/sqanti3_qc
mkdir -p $RO_OUTPUT_DIR
mkdir -p $FT_vs_RGC_OUTPUT_DIR

source activate /users/sparthib/.conda/envs/SQANTI3
SQANTI_DIR=~/SQANTI3-5.2.1

python $SQANTI_DIR/sqanti3_qc.py  ${ROs_GTF} ${REFERENCE_GTF} ${REFERENCE_GENOME_FASTA} \
    --skipORF -o ROs -d $RO_OUTPUT_DIR/sqanti3_qc --saturation \
    -t $SLURM_CPUS_PER_TASK --report skip --isoform_hits 
     ##positional arguments
     
python $SQANTI_DIR/sqanti3_qc.py  ${FT_vs_RGC_GTF} ${REFERENCE_GTF} ${REFERENCE_GENOME_FASTA} \
    --skipORF -o FT_vs_RGC -d $FT_vs_RGC_OUTPUT_DIR/sqanti3_qc --saturation \
    -t $SLURM_CPUS_PER_TASK --report skip --isoform_hits


echo "**** Job ends ****"
date +"%Y-%m-%d %T"
    