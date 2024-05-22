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

sample_wise_bed_dir=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output

#H9 FT vs HRGC
for file in $sample_wise_bed_dir/H9-FT_1_all_corrected.bed $sample_wise_bed_dir/H9-FT_2_all_corrected.bed $sample_wise_bed_dir/H9-hRGC_2_all_corrected.bed $sample_wise_bed_dir/H9-hRGC_1_all_corrected.bed;
do 
	cat "$file" >> $sample_wise_bed_dir/all_H9_FT_RGC_samples_corrected.bed
done

input=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output/all_H9_FT_RGC_samples_corrected.bed
output_dir=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output/all_H9_FT_RGC_sample_chr_split
mkdir -p $output_dir

for chr in {1..22} X Y M;
do
	echo $chr
	grep -w $chr $input > $output_dir/$chr.bed
done



#for files that have RO in their name (i.e. RO_1.bed, RO_2.bed, etc.)
for file in $sample_wise_bed_dir/*RO*.bed; 
do
	cat "$file" >> $sample_wise_bed_dir/all_RO_samples_corrected.bed
done

input=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output/all_RO_samples_corrected.bed
output_dir=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output/all_RO_sample_chr_split
mkdir -p $output_dir

for chr in {1..22} X Y M;
do
	echo $chr
	grep -w $chr $input > $output_dir/$chr.bed
done








echo "**** Job ends ****"
date +"%Y-%m-%d %T"



