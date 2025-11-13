#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 10
#SBATCH --job-name=download_VCFs
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/download_VCFs.txt
#SBATCH -e logs/download_VCFs.txt
#SBATCH -t 7-00:00:00



### download reference VCF files for BQRS

#similar to what Jianing did here 
#https://github.com/JianingYao/scRNA-seq-genetic-ancestry/blob/main/Preprocessing/GATK_pipeline/part1_bams.sh

known_sites=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/known_sites
mkdir -p $known_sites
cd $known_sites

wget -N https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -N https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

wget -N https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -N https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

wget -N https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -N https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi



