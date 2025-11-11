#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 15
#SBATCH --job-name=sqanti
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/sqanti.txt
#SBATCH -e logs/sqanti.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"

ENV_FILE="../../../.env"
if [ -f $ENV_FILE ]; then
    set -a
    source $ENV_FILE
    set +a
fi


#reference genome 
REFERENCE_GENOME_FASTA=$references_dir/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
REFERENCE_GTF=$references_dir/genome/GENCODE/primary_assembly/gencode.v46.basic.annotation.gtf    
#reference transcriptome 
REFERENCE_TRANSCRIPTOME=$references_dir/genome/GENCODE/primary_assembly/release_46_all_transcripts.fa

#sample transcriptome
BAMBU_GTF=$retina_lrs_dir/06_quantification/bambu/all_samples_extended_annotation_track_reads/extended_annotations_fa_contigs_only.gtf

# GTF (default): by default, SQANTI3 expects the transcriptome to be provided as a GTF file, 
# and we recommend to stick to this format if your transcriptome construction pipeline allows it.

OUTPUT_DIR=$retina_lrs_dir/06_quantification/bambu/all_samples_extended_annotation_track_reads/

source activate /users/sparthib/.conda/envs/SQANTI3
SQANTI_DIR=~/SQANTI3-5.2.1

python $SQANTI_DIR/sqanti3_qc.py  ${BAMBU_GTF} ${REFERENCE_GTF} ${REFERENCE_GENOME_FASTA} \
     -o all_samples -d $OUTPUT_DIR/sqanti3_qc --saturation \
    -t $SLURM_CPUS_PER_TASK --report skip --isoform_hits \
    --CAGE_peak $retina_lrs_dir/raw_data/cage/human.refTSS_v3.1.hg38.bed \
    --polyA_motif $retina_lrs_dir/raw_data/polya/mouse_and_human.polyA_motif.txt \
    --polyA_peak $retina_lrs_dir/raw_data/polya/atlas.clusters.2.0.GRCh38.96.bed 
     ##positional arguments


echo "**** Job ends ****"
date +"%Y-%m-%d %T"
    