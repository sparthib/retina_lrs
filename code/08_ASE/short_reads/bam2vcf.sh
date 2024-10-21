#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=bam2vcf.txt
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bam2vcf.txt
#SBATCH -e logs/bam2vcf.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"



ref_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa

# ml load samtools 
# samtools faidx $ref_fa

ml load gatk
# gatk CreateSequenceDictionary -R $ref_fa
# # https://samtools.github.io/bcftools/howtos/variant-calling.html

bam_files=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/filtered_bams/
output_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf

# bcftools mpileup -Ou --threads $SLURM_CPUS_PER_TASK \
# -f /dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa \
# -b /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/filtered_bams/bam_files.txt | bcftools call -mv -Ob > $output_dir/multi_sample.vcf

# gatk commands from here https://www.biostars.org/p/405702/
sample_names=(SRR1091088 SRR1091091)
# for sample in ${sample_names[@]}; do
#     echo "Processing $sample"
#     gatk --java-options "-Xmx4g" HaplotypeCaller \
#     -R $ref_fa \
#     -I $bam_files/$sample.sorted.bam \
#     -O $output_dir/$sample.g.vcf.gz \
#     -ERC GVCF
# done

echo "Combining GVCFs"
gatk --java-options "-Xmx96g -Xms96g" CombineGVCFs \
-R $reference \
--variant $output_dir/SRR1091088.g.vcf.gz \
--variant $output_dir/SRR1091091.g.vcf.gz \
-O $output_dir/combined.g.vcf.gz

echo "Genotyping"
gatk --java-options "-Xmx96g -Xms96g" GenotypeGVCFs \
-R $reference \
-V $output_dir/combined.g.vcf.gz \
-O $output_dir/combined.vcf.gz


echo "**** Job ends ****"
date +"%Y-%m-%d %T"




