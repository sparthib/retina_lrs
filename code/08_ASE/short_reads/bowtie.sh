#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=bowtie
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bowtie.txt
#SBATCH -e logs/bowtie.txt
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"

REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa.gz
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/sams
ml load bowtie/2.5.1

# Aligning paired reads
samples=(SRR1091088 SRR1091091 SRR1091092)

for sample in ${samples[@]}
do
    bowtie2 -x $REFERENCE_FASTA -p ${SLURM_CPUS_PER_TASK} \
    -1 example/reads/${sample}_1.fq -2 example/reads/${sample}_2.fq -S $OUTPUT_DIR/${sample}.sam
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

