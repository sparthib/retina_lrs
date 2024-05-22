#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 10
#SBATCH --job-name=flair_collapse
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flair_collapse_RO.%a.txt
#SBATCH -e logs/flair_collapse_RO.%a.txt
#SBATCH --array=1-25
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



CONFIG=/users/sparthib/retina_lrs/code/chromosome_number.txt
chr=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "chromosome" $chr

REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
REFERENCE_FASTQ=/dcs04/hicks/data/sparthib/retina_lrs/03_processed_fastqs/all_RO_files.fastq.gz
# longshot_output=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/longshot_vcfs/${sample}

#REFERENCE_FASTQ= list all files in the directory and get the file names
#concatenate all fq.gz files with RO in the name to all_fq_files.txt

source activate flair 

### take the bed files created from correction step 
bed_file=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output/all_RO_sample_chr_split/$chr.bed 

##output dir
collapsed_output=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/collapsed_output/all_RO_samples
mkdir -p $collapsed_output

flair collapse -r $REFERENCE_FASTQ -q $bed_file -g $REFERENCE_FASTA --threads 10 \
--stringent --check_splice --generate_map -f $REFERENCE_GTF --annotation_reliant generate \
-o $collapsed_output 
# --longshot_vcf $longshot_output/${sample}.vcf --longshot_bam $longshot_output/${sample}.bam \

echo "**** Job ends ****"
date +"%Y-%m-%d %T"