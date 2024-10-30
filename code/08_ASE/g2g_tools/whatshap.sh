#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --job-name=whatshap.txt
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/whatshap.txt
#SBATCH -e logs/whatshap.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ref_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
vcf_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf
source activate whatshap-env

whatshap phase -o $vcf_dir/SRR1091088.vcf --reference=$ref_fa $vcf_dir/SRR1091088.genotyped.vcf

echo "**** Job ends ****"
date +"%Y-%m-%d %T"



