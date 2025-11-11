#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=isoquant_FT_RGC
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/isoquant_FT_RGC.txt
#SBATCH -e logs/isoquant_FT_RGC.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

source activate isoquant 

ENV_FILE="../../.env"
if [ -f $ENV_FILE ]; then
    set -a
    source $ENV_FILE
    set +a
fi

sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "${sample}"
BAM_FOLDER=$retina_lrs_dir/05_bams/genome/primary_assembly/high_quality
REFERENCE_GTF=$references_dir/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf.gz
REFERENCE_FASTA=$references_dir/genome/GENCODE/primary_assembly/release_46_primary_genome.fa.gz

H9_FT_1=$BAM_FOLDER/H9-FT_1_primary_over_30_sorted.bam
H9_FT_2=$BAM_FOLDER/H9-FT_2_primary_over_30_sorted.bam
H9_hRGC_1=$BAM_FOLDER/H9-hRGC_1_primary_over_30_sorted.bam
H9_hRGC_2=$BAM_FOLDER/H9-hRGC_2_primary_over_30_sorted.bam

OUTPUT_FOLDER=$retina_lrs_dir/06_quantification/isoquant/high_quality/FT_RGC
mkdir -p $OUTPUT_FOLDER

isoquant.py --reference $REFERENCE_FASTA --data_type ont --genedb $REFERENCE_GTF \
--bam $H9_FT_1 $H9_FT_2 $H9_hRGC_1 $H9_hRGC_2  \
  --output $OUTPUT_FOLDER --clean_start --count_exons -t ${SLURM_CPUS_PER_TASK} \
  --complete_genedb --sqanti_output --check_canonical
  
  
conda deactivate

echo "**** Job ends ****"
date +"%Y-%m-%d %T"