#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --job-name=pfam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/pfam_bambu_FT_RGC.txt
#SBATCH -e logs/pfam_bambu_FT_RGC.txt
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"

source activate pfam 
cd /dcs04/hicks/data/sparthib/retina_lrs/PfamScan/pfam_scan
input_dir=/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC
./pfam_scan.py $input_dir/isoformSwitchAnalyzeR_isoform_AA_complete.fasta \
../ -out $input_dir/pfam_results.csv -cpu 20


conda deactivate 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

