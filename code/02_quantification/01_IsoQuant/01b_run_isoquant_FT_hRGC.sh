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

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "${sample}"
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa.gz

H9_FT_1=$BAM_FOLDER/H9-FT_1_primary_over_30_chr_only_sorted.bam
H9_FT_2=$BAM_FOLDER/H9-FT_2_primary_over_30_chr_only_sorted.bam
H9_hRGC_1=$BAM_FOLDER/H9-hRGC_1_primary_over_30_chr_only_sorted.bam
H9_hRGC_2=$BAM_FOLDER/H9-hRGC_2_primary_over_30_chr_only_sorted.bam


OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC
mkdir -p $OUTPUT_FOLDER

isoquant.py --reference $REFERENCE_FASTA --data_type ont --genedb $REFERENCE_GTF \
--bam $H9_FT_1 $H9_FT_2 $H9_hRGC_1 $H9_hRGC_2  \
  --output $OUTPUT_FOLDER --clean_start --count_exons -t ${SLURM_CPUS_PER_TASK} \
  --complete_genedb --sqanti_output --check_canonical
  
  
conda deactivate

echo "**** Job ends ****"
date +"%Y-%m-%d %T"