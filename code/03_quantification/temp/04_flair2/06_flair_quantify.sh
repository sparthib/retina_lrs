#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=flair_quantify
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flair_quantify.%a.txt
#SBATCH -e logs/flair_quantify.%a.txt
#SBATCH --array=1-15
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

READS_MANIFEST=/users/sparthib/retina_lrs/code/02_quantification/04_flair2/reads_manifest.tsv 
OUTDIR=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/quantify_output
ISOFORMS_FA=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/collapsed_output/$sample.isoforms.fa
ISOFORMS_BED=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/collapsed_output/$sample.isoforms.bed

source activate flair 

mkdir -p $OUTDIR/$sample

flair quantify -r $READS_MANIFEST -i $ISOFORMS_FA -o $OUTDIR/$sample --threads $SLURM_CPUS_PER_TASK \
--quality 30 --generate_map --check_splice --isoform_bed $ISOFORMS_BED

conda deactivate 


echo "**** Job ends ****"
date +"%Y-%m-%d %T"


