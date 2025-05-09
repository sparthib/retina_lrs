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
#SBATCH --array=1-8
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

ref_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
vcf_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/gvcf_ref_46
genome_bam_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality
whatshap_output_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output

source activate whatshap-env

# whatshap stats --gtf=$whatshap_output_dir/phased.gtf $whatshap_output_dir/phased.vcf

## input vcf file needs to be indexed prior to running whatshap
# used bgzip from htslib module for zipping vcf and then tabix to index
# ml load htslib
# bgzip -c $whatshap_output_dir/phased.vcf > $whatshap_output_dir/phased.vcf.gz
# tabix -p vcf $whatshap_output_dir/phased.vcf.gz

samples=(H9-BRN3B_hRO_2 H9-BRN3B-RO H9-CRX_hRO_2 H9-CRX_ROs_D45 H9-FT_1 H9-FT_2 H9-hRGC_1 H9-hRGC_2) 

lr_sample=${samples[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: $lr_sample"

whatshap haplotag -o $whatshap_output_dir/${lr_sample}.bam \
--reference $ref_fa $whatshap_output_dir/phased.vcf.gz $genome_bam_dir/${lr_sample}_primary_over_30_chr_only_sorted.bam \
--output-threads=19 --ignore-read-groups --sample=SRR1091091 --output-haplotag-list $whatshap_output_dir/${lr_sample}_haplotypes.tsv

whatshap split --output-h1 $whatshap_output_dir/${lr_sample}_h1.bam \
--output-h2 $whatshap_output_dir/${lr_sample}_h2.bam $whatshap_output_dir/${lr_sample}.bam $whatshap_output_dir/${lr_sample}_haplotypes.tsv



echo "**** Job ends ****"
date +"%Y-%m-%d %T"