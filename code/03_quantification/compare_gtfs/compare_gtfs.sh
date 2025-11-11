#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 4
#SBATCH --job-name=gffcompare
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --output=logs/gffcompare.out  # Standard output log
#SBATCH --error=logs/gffcompare.out   # Standard error log
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ENV_FILE="../../.env"
if [ -f $ENV_FILE ]; then
    set -a
    source $ENV_FILE
    set +a
fi

isoquant_gtf=$retina_lrs_code/06_quantification/isoquant/high_quality/all_samples/OUT/OUT.extended_annotation.gtf
bambu_gtf=$retina_lrs_code/06_quantification/bambu/all_samples_extended_annotation_track_reads/extended_annotations.gtf
gencode_gtf=$references_dir/genome/GENCODE/primary_assembly/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf
OUTPUT_DIR=$retina_lrs_code/06_quantification/gffcompare

mkdir -p $OUTPUT_DIR
~/gffcompare/gffcompare -R -r $gencode_gtf -o $OUTPUT_DIR/gencode_isoquant $isoquant_gtf
~/gffcompare/gffcompare -R -r $gencode_gtf -o $OUTPUT_DIR/gencode_bambu $bambu_gtf


~/gffcompare/gffcompare -R -r $bambu_gtf -o $OUTPUT_DIR/bambu_ref_isoquant_ext $isoquant_gtf

#only keep one instance of the same intron-chain between the two methods
~/gffcompare/gffcompare -R -r $bambu_gtf -o $OUTPUT_DIR/bambu_ref_gencode_ext_nodups $gencode_gtf -D
echo "**** Job ends ****"
date +"%Y-%m-%d %T"