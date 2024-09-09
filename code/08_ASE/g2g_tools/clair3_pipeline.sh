#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=clair3
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/clair.%a.txt
#SBATCH -e logs/clair.%a.txt
#SBATCH --array=9
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
stem_cell=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $5}' $CONFIG)
sample=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $2}' $CONFIG)

OUTPUT_DIR="/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/clair_vcfs/$stem_cell"
mkdir $OUTPUT_DIR
  

INPUT_DIR="/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only"        
MODEL_NAME="r941_prom_sup_g5014"
REFERENCE_FASTA="/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"

# source activate clair3
# 
# run_clair3.sh \
#   --bam_fn=$INPUT_DIR/${sample}_primary_over_30_chr_only_sorted.bam \                 ## change your bam file name here
#   --ref_fn=$REFERENCE_FASTA \                    ## change your reference file name here
#   --threads=20 \               ## maximum threads to be used
#   --platform="ont" \                   ## options: {ont,hifi,ilmn}
#   --model_path="/users/sparthib/.conda/envs/clair3/bin/models/${MODEL_NAME}" \ 
#   --output=${OUTPUT_DIR}     ## output path prefix 

ml load singularity

singularity pull docker://hkubal/clair3:latest

# run clair3 like this afterward
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  clair3_latest.sif \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/${sample}_primary_over_30_chr_only_sorted.bam \
  --ref_fn=${REFERENCE_FASTA} \
  --threads=20 \
  --platform="ont" \
  --model_path="/users/sparthib/.conda/envs/clair3/bin/models/r941_prom_sup_g5014" \
  --output=${OUTPUT_DIR}
  
  
echo "**** Job ends ****"
date +"%Y-%m-%d %T"
