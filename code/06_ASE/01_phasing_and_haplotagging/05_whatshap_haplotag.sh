#!/bin/bash

#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=haplotag
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/haplotag.%a.txt
#SBATCH -e logs/haplotag.%a.txt
#SBATCH --array=1-11
#SBATCH --time=7-00:00:00

### whatshap phases the variants we found using GATK with the help 
### of our long-read BAMs
### This phased VCF is further used to haplotag our BAM files. 


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

ref_fa=$references_dirgenome/GENCODE/primary_assembly/release_46_primary_genome.fa
vcf_dir=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/gvcf_ref_46
genome_bam_dir=$retina_lrs_dir/05_bams/genome/primary_assembly/high_quality
phased_vcf=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/whatshap_output/all_samples_phased.vcf.gz
phased_vcf_H9_EP1=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/whatshap_output/all_samples_H9_and_EP1_phased.vcf.gz
whatshap_output_dir=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/whatshap_output_single_sample
whatshap_output_dir_H9_EP1=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1
#mkdir -p $whatshap_output_dir_H9_EP1

source activate whatshap-env

# whatshap stats --gtf=$whatshap_output_dir/phased.gtf $whatshap_output_dir/phased.vcf

## input vcf file needs to be indexed prior to running whatshap
## used bgzip from htslib module for zipping vcf and then tabix to index
#ml load htslib
# bgzip -c $phased_vcf_H9_EP1 > $phased_vcf_H9_EP1.gz
# tabix -p vcf $phased_vcf_H9_EP1.gz

samples=(H9-BRN3B_hRO_2 H9-BRN3B-RO H9-CRX_hRO_2 H9-CRX_ROs_D45 H9-FT_1 H9-FT_2 H9-hRGC_1 H9-hRGC_2 EP1-BRN3B-RO EP1-WT_hRO_2 EP1-WT_ROs_D45) 
# Stage 2, Stage 3, Stage 2, Stage1, FT1, FT2, RGC1, RGC2, Stage 3, Stage 2, Stage 1 
lr_sample=${samples[$SLURM_ARRAY_TASK_ID - 1]}

echo "**** Haplotagging sample: $lr_sample ****"
whatshap haplotag -o $whatshap_output_dir_H9_EP1/${lr_sample}.bam \
--reference $ref_fa $phased_vcf_H9_EP1 $genome_bam_dir/${lr_sample}_primary_over_30_chr_only_sorted.bam \
--output-threads=19 --ignore-read-groups --output-haplotag-list $whatshap_output_dir_H9_EP1/${lr_sample}_haplotypes.tsv

echo "**** Splitting haplotagged BAM into haplotypes ****"
whatshap split --output-h1 $whatshap_output_dir_H9_EP1/${lr_sample}_h1.bam \
--output-h2 $whatshap_output_dir_H9_EP1/${lr_sample}_h2.bam $whatshap_output_dir_H9_EP1/${lr_sample}.bam $whatshap_output_dir_H9_EP1/${lr_sample}_haplotypes.tsv


echo "**** Job ends ****"
date +"%Y-%m-%d %T"