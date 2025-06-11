#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 10
#SBATCH --job-name=igv
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/merge_bams.txt
#SBATCH -e logs/merge_bams.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

output_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/igv

ml load samtools 


for dir in SYNCRIP BAK1 BSG; do
   samtools merge $output_dir/${dir}/Stage_3_${dir}.bam $output_dir/${dir}/EP1-BRN3B-RO_$dir.bam $output_dir/$dir/H9-BRN3B-RO_$dir.bam
   samtools sort -@ 10 -o $output_dir/${dir}/Stage_3_${dir}_sorted.bam $output_dir/${dir}/Stage_3_${dir}.bam
   samtools index $output_dir/${dir}/Stage_3_${dir}_sorted.bam > $output_dir/${dir}/Stage_3_${dir}.bam.bai
   
   samtools merge $output_dir/${dir}/Stage_2_${dir}.bam $output_dir/${dir}/EP1-WT_hRO_2_$dir.bam $output_dir/$dir/H9-BRN3B_hRO_2_$dir.bam $output_dir/$dir/H9-CRX_hRO_2_$dir.bam
   samtools sort -@ 10 -o $output_dir/${dir}/Stage_2_${dir}_sorted.bam $output_dir/${dir}/Stage_2_${dir}.bam
   samtools index $output_dir/${dir}/Stage_2_${dir}_sorted.bam > $output_dir/${dir}/Stage_2_${dir}.bam.bai
   
   samtools merge $output_dir/${dir}/Stage_1_${dir}.bam $output_dir/${dir}/EP1-WT_ROs_D45_$dir.bam $output_dir/$dir/H9-CRX_ROs_D45_$dir.bam 
   samtools sort -@ 10 -o $output_dir/${dir}/Stage_1_${dir}_sorted.bam $output_dir/${dir}/Stage_1_${dir}.bam
   samtools index $output_dir/${dir}/Stage_1_${dir}_sorted.bam > $output_dir/${dir}/Stage_1_${dir}.bam.bai
   
   samtools merge $output_dir/${dir}/FT_${dir}.bam $output_dir/${dir}/H9-FT_1_$dir.bam $output_dir/$dir/H9-FT_2_$dir.bam
   samtools sort -@ 10 -o $output_dir/${dir}/FT_${dir}_sorted.bam $output_dir/${dir}/FT_${dir}.bam
   samtools index $output_dir/${dir}/FT_${dir}_sorted.bam > $output_dir/${dir}/FT_${dir}.bam.bai
   
   samtools merge $output_dir/${dir}/RGC_${dir}.bam $output_dir/${dir}/H9-hRGC_1_$dir.bam $output_dir/$dir/H9-hRGC_2_$dir.bam
   samtools sort -@ 10 -o $output_dir/${dir}/RGC_${dir}_sorted.bam $output_dir/${dir}/RGC_${dir}.bam
   samtools index $output_dir/${dir}/RGC_${dir}_sorted.bam > $output_dir/${dir}/RGC_${dir}.bam.bai
   
done