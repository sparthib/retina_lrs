#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=flair_collapse
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flair_collapse.%a.txt
#SBATCH -e logs/flair_collapse.%a.txt
#SBATCH --array=1-15
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"



CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
REFERENCE_FASTQ=/dcs04/hicks/data/sparthib/retina_lrs/03_processed_fastqs/${sample}.fastq.gz
longshot_output=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/longshot_vcfs/${sample}


source activate flair 

### take the bed files created from correction step 
bed_file=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output/${sample}_all_corrected.bed
collapsed_output=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/collapsed_output/${sample}

mkdir -p $collapsed_output

flair collapse -r $REFERENCE_FASTQ -q $bed_file -g $REFERENCE_FASTA --threads 20 \
-stringent --check_splice --generate_map --annotation_reliant generate \
# --longshot_vcf $longshot_output/${sample}.vcf --longshot_bam $longshot_output/${sample}.bam \
-o $collapsed_output 


echo "**** Job ends ****"
date +"%Y-%m-%d %T"