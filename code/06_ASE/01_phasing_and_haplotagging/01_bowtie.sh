#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=200G
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

ENV_FILE="../../../.env"
if [ -f $ENV_FILE ]; then
    set -a
    source $ENV_FILE
    set +a
fi

INPUT_DIR=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/fastq
OUTPUT_DIR=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/sams_ref_46
mkdir -p $OUTPUT_DIR
##../H9_DNA_Seq_data/sams has outputs from aligning to GRCh38_noalt_as


ml load bowtie/2.5.1

ref_fa=$references_dir/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
ref_fa_index_dir=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/release_46_index
mkdir -p $ref_fa_index_dir
bowtie2-build $ref_fa $ref_fa_index_dir

#BT_INDEX_PATH=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/GRCh38_noalt_as/
#BT_INDEX_PATH was downloaded from bowtie2 reference

# Aligning paired reads
samples=(SRR1091088 SRR1091091 SRR1091092)

for sample in ${samples[@]}
do
    bowtie2 -x $ref_fa_index_dir -p ${SLURM_CPUS_PER_TASK} \
    -1 $INPUT_DIR/${sample}_1.fastq -2 $INPUT_DIR/${sample}_2.fastq -S $OUTPUT_DIR/${sample}.sam \
    --rg-id ${sample} --rg "SM:${sample}" --rg "PL:ILLUMINA" 
    echo "**** done aligning $sample ****"
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

