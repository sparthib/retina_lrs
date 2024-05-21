#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=split_bed
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/split_bed.txt
#SBATCH -e logs/split_bed.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"



input=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output/all_samples_corrected.bed

for chr in `cut -f 1 $input | sort | uniq`;
do
	echo $chr
	grep -w $chr $input > $chr.bed
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"