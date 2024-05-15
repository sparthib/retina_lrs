#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=flair_assign_var
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flair_assign_variant.%a.txt
#SBATCH -e logs/flair_assign_variant.%a.txt
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

REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
REFERENCE_FASTQ=/dcs04/hicks/data/sparthib/retina_lrs/03_processed_fastqs/${sample}.fastq.gz
longshot_output=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/longshot_vcfs/${sample}/${sample}.vcf
genome_bam=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/${sample}_sorted.bam
out_bed=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/assign_variants/${sample}.bed
out_map=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/assign_variants/${sample}.map
out_vcf=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/assign_variants/${sample}.vcf
source activate flair 

### take the bed files created from correction step 
bed_file=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output/${sample}_all_corrected.bed
collapsed_output=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/collapsed_output/${sample}

mkdir -p $collapsed_output

assign_variants_to_transcripts –bam $genome_bam -i $sample.isoforms.bed \
-v $longshot_output –map $sample.isoform.read.map.txt \
–bed_out $out_bed –map_out $out_map > $out_vcf


echo "**** Job ends ****"
date +"%Y-%m-%d %T"