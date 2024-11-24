#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 4
#SBATCH --job-name=gffcompare
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --output=logs/gffcompare.out  # Standard output log
#SBATCH --error=logs/gffcompare.err   # Standard error log
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


isoquant_gtf=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT/OUT.extended_annotation.gtf
bambu_gtf=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads/extended_annotations.gtf
gencode_gtf=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/gffcompare

mkdir -p $OUTPUT_DIR
/users/sparthib/gffcompare/gffcompare -R -r $gencode_gtf -o $OUTPUT_DIR/gencode_isoquant $isoquant_gtf
/users/sparthib/gffcompare/gffcompare -R -r $gencode_gtf -o $OUTPUT_DIR/gencode_bambu $bambu_gtf

echo "**** Job ends ****"
date +"%Y-%m-%d %T"