#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=sqanti
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/sqanti.%a.txt
#SBATCH -e logs/sqanti.%a.txt
#SBATCH --array=1
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"



CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

#reference genome 
REFERENCE_GENOME_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
#reference transcriptome 
REFERENCE_TRANSCRIPTOME=/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts_short_header.fa
#sample transcriptome
SAMPLE_TRANSCRIPTOME=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/collapsed_output/${sample}.isoforms.fa
SAMPLE_GTF=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/collapsed_output/${sample}.isoforms.gtf
# GTF (default): by default, SQANTI3 expects the transcriptome to be provided as a GTF file, 
# and we recommend to stick to this format if your transcriptome construction pipeline allows it
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/sqanti_output

source activate /users/sparthib/.conda/envs/SQANTI3
SQANTI_DIR=~/SQANTI3-5.2.1

python $SQANTI_DIR/sqanti3_qc.py   \
    --skipORF --fasta ${SAMPLE_TRANSCRIPTOME} \
    -o $sample --dir $OUTPUT_DIR ${sample}_sqanti3_qc --saturation \
    -t $SLURM_CPUS_PER_TASK --report pdf --isoform_hits \
    ${SAMPLE_GTF} ${REFERENCE_GTF} ${REFERENCE_GENOME_FASTA} ##positional arguments


echo "**** Job ends ****"
date +"%Y-%m-%d %T"
    