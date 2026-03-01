#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 10
#SBATCH --job-name=subset_hap_bams
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-11
#SBATCH -o logs/subset_hap_bams.%a.txt
#SBATCH -e logs/subset_hap_bams.%a.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID}"


# Samples match order used in 05_whatshap_haplotag.sh
samples=(H9-BRN3B_hRO_2 H9-BRN3B-RO H9-CRX_hRO_2 H9-CRX_ROs_D45 \
         H9-FT_1 H9-FT_2 H9-hRGC_1 H9-hRGC_2 \
         EP1-BRN3B-RO EP1-WT_hRO_2 EP1-WT_ROs_D45)
sample=${samples[$SLURM_ARRAY_TASK_ID - 1]}
echo "Processing sample: $sample"

bam_dir=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1
output_dir=$bam_dir/igv
bed_file=$output_dir/genes_of_interest.bed   # written by 06_genes_of_interest_to_bed.R

mkdir -p $output_dir

ml load samtools

for hap in h1 h2; do
    input_bam=$bam_dir/${sample}_${hap}.bam

    # Sort and index if index is missing
    if [ ! -f ${input_bam}.bai ]; then
        echo "Sorting and indexing $input_bam"
        samtools sort -@ 10 -o ${input_bam%.bam}_sorted_tmp.bam $input_bam
        mv ${input_bam%.bam}_sorted_tmp.bam $input_bam
        samtools index $input_bam
    fi

    # Subset to each gene region listed in the BED file
    while IFS=$'\t' read -r chrom start end gene_name rest; do
        mkdir -p $output_dir/$gene_name

        # BED is 0-based; samtools view region is 1-based inclusive
        region="${chrom}:$((start + 1))-${end}"
        out_bam=$output_dir/$gene_name/${sample}_${hap}_${gene_name}.bam

        echo "  [$hap] $gene_name ($region)"
        samtools view -b $input_bam "$region" > $out_bam
        samtools index $out_bam
    done < $bed_file

done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
