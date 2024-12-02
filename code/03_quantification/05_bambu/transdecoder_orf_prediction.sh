#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 15
#SBATCH --job-name=transdecoder
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/transdecoder.txt
#SBATCH -e logs/transdecoder.txt
#SBATCH --time=7-00:00:00


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "****"

BAMBU_DIR=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads
REFERENCE_GENOME_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
BAMBU_GTF=$BAMBU_DIR/extended_annotations.gtf
# cut -f1 $BAMBU_GTF | sort | uniq > ~/gtf_contigs.txt
# grep "^>" $REFERENCE_GENOME_FASTA | cut -d' ' -f1 | sed 's/>//' > ~/fasta_contigs.txt
# sort ~/fasta_contigs.txt -o ~/fasta_contigs.txt
# comm -23 ~/gtf_contigs.txt ~/fasta_contigs.txt > ~/missing_contigs.txt
# sort ~/missing_contigs.txt -o ~/missing_contigs_sorted.txt

grep -F -v -f ~/missing_contigs_sorted.txt $BAMBU_GTF > $BAMBU_DIR/extended_annotations_fa_contigs_only.gtf

cd /users/sparthib/TransDecoder-TransDecoder-v5.7.1

./util/gtf_genome_to_cdna_fasta.pl $BAMBU_DIR/extended_annotations_fa_contigs_only.gtf $REFERENCE_GENOME_FASTA > $BAMBU_DIR/extended_annotations.fasta

./util/gtf_to_alignment_gff3.pl $BAMBU_DIR/extended_annotations_fa_contigs_only.gtf > $BAMBU_DIR/extended_annotations_fa_contigs_only.gtf.gff3


./TransDecoder.LongOrfs -t $BAMBU_DIR/extended_annotations.fasta

util/cdna_alignment_orf_to_genome_orf.pl $BAMBU_DIR/extended





