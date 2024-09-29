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
#SBATCH --array=1-3


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID}"
echo "****"



#reference genome 
REFERENCE_GENOME_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
#reference transcriptome 
REFERENCE_TRANSCRIPTOME=/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts_short_header.fa

samples=("10x_D100-EP1" "10x_D200-EP1-1" "10x_D200-EP1-2")
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

#sample transcriptome
GTF=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/$sample/OUT/OUT.extended_annotation.gtf 

OUTPUT_DIR=//dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/$sample/OUT/sqanti3_qc
mkdir -p $OUTPUT_DIR

source activate /users/sparthib/.conda/envs/SQANTI3
SQANTI_DIR=~/SQANTI3-5.2.1

python $SQANTI_DIR/sqanti3_qc.py $GTF $REFERENCE_GTF $REFERENCE_GENOME_FASTA \
     -o ROs -d $RO_OUTPUT_DIR/sqanti3_qc --saturation \
    -t $SLURM_CPUS_PER_TASK --report skip --isoform_hits \
    --CAGE_peak /users/sparthib/retina_lrs/raw_data/cage/human.refTSS_v3.1.hg38.bed \
    --polyA_motif /users/sparthib/retina_lrs/raw_data/polya/mouse_and_human.polyA_motif.txt \
    --polyA_peak /users/sparthib/retina_lrs/raw_data/polya/atlas.clusters.2.0.GRCh38.96.bed 
     

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
    