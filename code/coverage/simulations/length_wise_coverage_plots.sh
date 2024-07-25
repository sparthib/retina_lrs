
transcripts=/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/FT_vs_RGC/long_transcripts.txt
input=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/supplementary_filtered/H9-FT_1_sorted.bam  


#extract reads that match ensemble transcript ids
samtools view -h $input | awk -F'\t' 'NR==FNR {isoforms[$2]; next} ($3 in isoforms) || ($1 ~ /^@/)' $transcripts - | \
samtools view -Sb -h > H9_FT_1_long_transcripts.bam


#get a list of reads from the transcriptome aligned bam files 
long_transcriptome_bam=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/supplementary_filtered/H9_FT_1_long_transcripts.bam
reads=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/supplementary_filtered/long_transcript_read_ids.txt 
samtools view -h $long_transcriptome_bam | awk -F'\t' '($1 ~ /^@/) {print $1} ($1 !~ /^@/) {print $1}' > $reads

genome_bam_input=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9-FT_1_primary_over_30_chr_only_sorted.bam


#extract reads based on read id
samtools view -h $genome_bam_input | awk -F'\t' 'NR==FNR {read[$1]; next} ($1 in read) || ($1 ~ /^@/)' $reads - | \
samtools view -Sb -h > /dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9_FT_1_long_transcripts.bam
#sort and index

samtools index /dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9_FT_1_short_transcripts.bam
samtools index /dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9_FT_1_medium_transcripts.bam
samtools index /dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9_FT_1_long_transcripts.bam 
