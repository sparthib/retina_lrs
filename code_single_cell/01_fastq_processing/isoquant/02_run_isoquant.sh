#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=isoquant
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/isoquant.%a.txt
#SBATCH -e logs/isoquant.%a.txt
#SBATCH --array=4
#SBATCH --time=7-00:00:00

#try running for all chromosomes

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


source activate isoquant 

CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo $sample


BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only

REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa.gz

OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/$sample
mkdir -p $OUTPUT_FOLDER

isoquant.py -d ont --bam ${BAM_FOLDER}/${sample}_with_tags.bam \
--reference $REFERENCE_FASTA --genedb $REFERENCE_GTF --complete_genedb \
--output $OUTPUT_FOLDER -t ${SLURM_CPUS_PER_TASK} \
--sqanti_output --check_canonical --count_exons --bam_tags CB \
--model_construction_strategy default_ont --report_canonical auto --read_group tag:CB

conda deactivate

echo "**** Job ends ****"
date +"%Y-%m-%d %T"


