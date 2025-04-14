#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=300G
#SBATCH -c 10
#SBATCH --job-name=scnanoseq
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/scnanoseq.txt
#SBATCH -e logs/scnanoseq.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


ml load nextflow/24.10.5
ml load singularity

ref_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
ref_gtf=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf
ref_transcript_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_all_transcripts_short_header.fa

output_dir=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/scnanoseq

nextflow run nf-core/scnanoseq \
   -profile singularity \
   --input /users/sparthib/retina_lrs/code_single_cell/01_fastq_processing/scnanoseq/samplesheet.csv \
   --outdir $output_dir \
   -work-dir $output_dir/work \
   --genome_fasta $ref_fa \
   --transcript_fasta $ref_transcript_fa \
   --gtf $ref_gtf \
   --barcode_format 10X_3v3 \
   --dedup_tool umi_tools \
   --quantifier isoquant,oarfish
   


echo "**** Job ends ****"
date +"%Y-%m-%d %T"