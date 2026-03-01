#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 10
#SBATCH --job-name=merge_hap_bams
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/merge_hap_bams.txt
#SBATCH -e logs/merge_hap_bams.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

retina_lrs_dir=/dcs04/hicks/data/sparthib/retina_lrs
output_dir=$retina_lrs_dir/09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1/igv
bed_file=$output_dir/genes_of_interest.bed   # written by 06_genes_of_interest_to_bed.R

ml load samtools

# Helper: merge, sort, index; then remove unsorted intermediate
merge_sort_index() {
    local out_merged=$1
    local out_sorted=${out_merged%.bam}_sorted.bam
    shift
    samtools merge -f $out_merged "$@"
    samtools sort -@ 10 -o $out_sorted $out_merged
    samtools index $out_sorted
    rm $out_merged
}

while IFS=$'\t' read -r chrom start end gene_name rest; do
    echo "Merging gene: $gene_name"
    d=$output_dir/$gene_name

    for hap in h1 h2; do
        HAP="${hap/h/H}"   # h1 -> H1, h2 -> H2

        # Stage 3 (BRN3B stage): EP1-BRN3B-RO + H9-BRN3B-RO
        merge_sort_index \
            $d/Stage_3_${HAP}_${gene_name}.bam \
            $d/EP1-BRN3B-RO_${hap}_${gene_name}.bam \
            $d/H9-BRN3B-RO_${hap}_${gene_name}.bam

        # Stage 2 (hRO_2): EP1-WT_hRO_2 + H9-BRN3B_hRO_2 + H9-CRX_hRO_2
        merge_sort_index \
            $d/Stage_2_${HAP}_${gene_name}.bam \
            $d/EP1-WT_hRO_2_${hap}_${gene_name}.bam \
            $d/H9-BRN3B_hRO_2_${hap}_${gene_name}.bam \
            $d/H9-CRX_hRO_2_${hap}_${gene_name}.bam

        # Stage 1 (D45): EP1-WT_ROs_D45 + H9-CRX_ROs_D45
        merge_sort_index \
            $d/Stage_1_${HAP}_${gene_name}.bam \
            $d/EP1-WT_ROs_D45_${hap}_${gene_name}.bam \
            $d/H9-CRX_ROs_D45_${hap}_${gene_name}.bam

        # FT: H9-FT_1 + H9-FT_2
        merge_sort_index \
            $d/FT_${HAP}_${gene_name}.bam \
            $d/H9-FT_1_${hap}_${gene_name}.bam \
            $d/H9-FT_2_${hap}_${gene_name}.bam

        # RGC: H9-hRGC_1 + H9-hRGC_2
        merge_sort_index \
            $d/RGC_${HAP}_${gene_name}.bam \
            $d/H9-hRGC_1_${hap}_${gene_name}.bam \
            $d/H9-hRGC_2_${hap}_${gene_name}.bam

    done

done < $bed_file

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
