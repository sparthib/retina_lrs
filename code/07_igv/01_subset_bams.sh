#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 10
#SBATCH --job-name=igv
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1,2,5-6,9-15
#SBATCH -o logs/igv_SYNCRIP.%a.txt
#SBATCH -e logs/igv_SYNCRIP.%a.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID}"


CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

bam_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality
output_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/igv

mkdir -p $output_dir

# Load the required module
ml load samtools 

#REEP6 chromosome 	NC_000019.10 (1491181..1497927)
#ADD1 chromosome 	NC_000004.12 (2843844..2930062)
#SYNCRIP chr6:85613976-85642991
#BAK1 chr6:33572547-33580293
#BSG  chr19:571277-583494
#"PROM1", "RP1","CRB1","CRX"
#RP1 8: 54,509,422-54,871,720
#CRB1 1: 197,268,204-197,478,455
#CRX

for dir in SYNCRIP BAK1 BSG PROM1 RP1 CRB1 CRX; do
    mkdir -p $output_dir/$dir
done

# samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr19:1491181-1497927" > $output_dir/REEP6/${sample}_REEP6.bam
# samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr4:2843844-2930062" > $output_dir/ADD1/${sample}_ADD1.bam
samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr6:85613976-85642991" > $output_dir/SYNCRIP/${sample}_SYNCRIP.bam
samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr6:33572547-33580293" > $output_dir/BAK1/${sample}_BAK1.bam
samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr19:571277-583494" > $output_dir/BSG/${sample}_BSG.bam
samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr4:2850000-2855000" > $output_dir/PROM1/${sample}_PROM1.bam
samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr8:54509422-54871720" > $output_dir/RP1/${sample}_RP1.bam
samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr1:197268204-197478455" > $output_dir/CRB1/${sample}_CRB1.bam
samtools view -b $bam_dir/${sample}_primary_over_30_chr_only_sorted.bam "chr19:47819779-47843330" > $output_dir/CRX/${sample}_CRX.bam

## index bam files
for dir in SYNCRIP BAK1 BSG PROM1 RP1 CRB1 CRX; do
    samtools index $output_dir/$dir/${sample}_$dir.bam > $output_dir/$dir/${sample}_$dir.bam.bai
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

