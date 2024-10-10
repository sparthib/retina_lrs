#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=filter_vcfs
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/filter_vcfs.txt
#SBATCH -e logs/filter_vcfs.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"




ml load gatk

output_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf

#params based on this tutorial https://sib-swiss.github.io/NGS-variants-training/2021.9/day2/filtering_evaluation/
#split VCFs into SNPs and Indels

gatk SelectVariants \
--variant $output_dir/multi_sample.vcf \
-select-type SNP \
--output $output_dir/SNP.vcf

gatk SelectVariants \
--variant $output_dir/multi_sample.vcf \
-select-type INDEL \
--output $output_dir/INDEL.vcf

gatk VariantFiltration \
--variant $output_dir/SNP.vcf\
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "SOR > 3.0" --filter-name "SOR3" \
--filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--output $output_dir/filtered_SNP.vcf


gatk VariantFiltration \
--variant $output_dir/INDEL.vcf \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "FS > 200.0" --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
--output $output_dir/filtered_INDEL.vcf


echo "**** Job ends ****"
date +"%Y-%m-%d %T"