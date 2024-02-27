

#######.LongQCs  ######

OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/03_processed_fastqs/H9-CRX_ROs_D45_LongQC/
SAMPLE=H9-CRX_ROs_D45.fastq.gz
INPUT_SAMPLE=/dcs04/hicks/data/sparthib/retina_lrs/03_processed_fastqs/$SAMPLE



H9-CRX_ROs_D45_LongQC/trimmed

python ~/LongQC/longQC.py --ont -o $OUTPUT_DIR $INPUT_SAMPLE \
--adapter_5 "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT" --adapter_3 "GCAATACGTAACTGAACGAAGTACAGG" \
-p 15 -d -c $OUTPUT_DIR/trimmed




######## Qualimap. ########
#qualimap rnaseq -bam YZ_15T_hRGC_sorted.chr21.bam -gtf /dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf -outformat PDF
  
##ALGORITHM 
### Interested genes list 

# suppose it is of the format: 

#gene_name, chromosome:  chr4 , GRange: chr4:325000-400000, sample_name, 
# load file name: that's something like "../$sample/$chr.bam 



######### bam2plot ###########
PATH_TO_BAM=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/hRGC_sorted.bam
OUTPUT=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/TEST
bam2plot --bam $PATH_TO_BAM --sample_name hRGC_bam2plot --outpath $OUTPUT  --threshold 0 -i False -s False

####### RSeQC ##########
# /users/sparthib/.conda/envs/rseqc/bin/geneBody_coverage.py
REF_BED=/dcs04/hicks/data/sparthib/references/rseqc_housekeeping.bed
BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/YZ_15T_hRGC_sorted.chr21.bam
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/rseqc_test
/users/sparthib/.conda/envs/rseqc/bin/geneBody_coverage.py -r $REF_BED -i $BAM_FILE -o $OUTPUT_DIR

######karyoploteR


for file in ./*.bam 
do 
  samtools sort $file -o ./sorted/${file}_sorted.bam
# rm ${BAM_FOLDER}/${sample}.bam
  samtools index ./sorted/${file}_sorted.bam ./sorted/${file}_sorted.bam.bai
done




######### QoRTS #########
JAR_PATH=/users/sparthib/QoRTS/QoRTs.jar
BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/YZ_15T_hRGC_sorted.chr21.bam
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/QoRTS_test/
GTF_ANNO=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz

java -Xmx35G -jar $JAR_PATH QC --generatePlots --singleEnded --minMAPQ 60 $BAM_FILE $GTF_ANNO $OUTPUT_DIR 
## all the file arguments have to be in the end. 


#### CuteSV #####
BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/YZ_15T_hRGC_sorted.chr21.bam
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/TEST
cuteSV -S YZ_15T_hRGC_sorted.chr21 --retain_work_dir --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 $BAM_FILE $REFERENCE_FASTA CUTESV_Test.vcf $OUTPUT_DIR


##### LongGF #####
cd ~/LongGF/bin
BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/YZ_15T_hRGC_sorted.chr21.bam
GTF_ANNO=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
OUTPUT_DIR=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/TEST
./LongGF $BAM_FILE $GTF_ANNO 100 50 200 > $OUTPUT_DIR/LongGF.run.on.FEB.13.log



#### JBrowse Export ####
FASTA=/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts_short_header.fa
BAM=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/MAPQ_filtered
OUT=/users/sparthib/retina_lrs/plots/coverage/JB2Export
jb2export --fasta $FASTA --bam $BAM/EP1-BRN3B-RO_sorted.bam snpcov force:true height:600 --loc ENST00000333511.9 --out $OUT/EP1-BRN3B-RO_sorted_BSG.svg 

#### BEDTOOLS BAMTOBED ####
### GO TO ~/bedtools2 first 
cd ~/bedtools2
BAM=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/MAPQ_FILTERED/EP1-BRN3B-RO_sorted.bam
OUTPUT=/dcs04/hicks/data/sparthib/retina_lrs/05b_beds
bedtools bamtobed -i $BAM > $OUTPUT/EP1-BRN3B-RO_sorted.bed


