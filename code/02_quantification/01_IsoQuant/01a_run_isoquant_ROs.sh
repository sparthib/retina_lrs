#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=200G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=isoquant_RO
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/isoquant_ROs.txt
#SBATCH -e logs/isoquant_ROs.txt
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

EP1_BRN3B_RO=$BAM_FOLDER/EP1-BRN3B-RO_primary_over_30_chr_only_sorted.bam
EP1_WT_hRO_2=$BAM_FOLDER/EP1-WT_hRO_2_primary_over_30_chr_only_sorted.bam
EP1_WT_ROs_D45=$BAM_FOLDER/EP1-WT_ROs_D45_primary_over_30_chr_only_sorted.bam
H9_BRN3B_hRO_2=$BAM_FOLDER/H9-BRN3B_hRO_2_primary_over_30_chr_only_sorted.bam
H9_BRN3B_RO=$BAM_FOLDER/H9-BRN3B-RO_primary_over_30_chr_only_sorted.bam
H9_CRX_hRO_2=$BAM_FOLDER/H9-CRX_hRO_2_primary_over_30_chr_only_sorted.bam
H9_CRX_ROs_D45=$BAM_FOLDER/H9-CRX_ROs_D45_primary_over_30_chr_only_sorted.bam

OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/ROs
mkdir -p $OUTPUT_FOLDER

isoquant.py --reference $REFERENCE_FASTA --data_type ont --genedb $REFERENCE_GTF \
--bam $EP1_BRN3B_RO $EP1_WT_hRO_2 $EP1_WT_ROs_D45 $H9_BRN3B_hRO_2 $H9_BRN3B_RO $H9_CRX_hRO_2 $H9_CRX_ROs_D45 \
  --output $OUTPUT_FOLDER --clean_start --count_exons -t ${SLURM_CPUS_PER_TASK} --complete_genedb \
  --sqanti_output --check_canonical 
  
  
conda deactivate

echo "**** Job ends ****"
date +"%Y-%m-%d %T"