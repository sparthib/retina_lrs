#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=350G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=isoquant
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/isoquant_all_samples
#SBATCH -e logs/isoquant_all_samples
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

source activate isoquant 

# samples=("10x_D100-EP1" "10x_D200-EP1-1" "10x_D200-EP1-2")
# sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only/tagged

mapfile -t files < <(find "$BAM_FOLDER" -name "*_deduped_sorted.bam")

rep1="${files[0]}"
rep2="${files[1]}"
rep3="${files[2]}"
rep4="${files[3]}"
rep5="${files[4]}"
rep6="${files[5]}"
rep7="${files[6]}"
rep8="${files[7]}"
rep9="${files[8]}"
rep10="${files[9]}"
rep11="${files[10]}"
rep12="${files[11]}"

REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa.gz

OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/all_samples
mkdir -p $OUTPUT_FOLDER

isoquant.py -d ont --bam $rep1 $rep2 $rep3 $rep4 $rep5 $rep6 $rep7 $rep8 $rep9 $rep10 $rep11 $rep12 \
--reference $REFERENCE_FASTA --genedb $REFERENCE_GTF --complete_genedb \
--output $OUTPUT_FOLDER -t ${SLURM_CPUS_PER_TASK} \
--sqanti_output --check_canonical --count_exons --bam_tags CB \
--model_construction_strategy default_ont --report_canonical auto --read_group tag:CB

conda deactivate

echo "**** Job ends ****"
date +"%Y-%m-%d %T"