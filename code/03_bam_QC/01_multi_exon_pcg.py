import pysam
from collections import Counter
import pandas as pd
import numpy as np
import os
import pyranges as pr
import sys

#from https://github.molgen.mpg.de/MayerGroup/NGN3_paper_code/blob/main/ONT_seq_analysis/long_short_comparison/junction_per_read.py
# Load the GTF file
gtf_fn = "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
gtf = pr.read_gtf(gtf_fn, as_df = True)

curr_gene = ''
goi = list()
curr_num_exon = 0
for i in range(gtf.shape[0]):
    if gtf.iloc[i,2] == 'gene' and gtf.iloc[i,9] == 'protein_coding':
        if curr_num_exon > 1:
            goi.append(curr_gene)
        curr_gene = gtf.iloc[i,8]
        curr_num_exon = 0
    if gtf.iloc[i,2] == 'exon' and gtf.iloc[i,8] == curr_gene:
        curr_num_exon += 1


genic_gtf = gtf[gtf['Feature'] == 'gene']
genic_gtf = genic_gtf[genic_gtf['gene_id'].isin(goi)]
genic_gtf = genic_gtf.iloc[:, 0:14]

###Get bam file from command line arg
if len(sys.argv) > 1:
    sample = sys.argv[1]
    print(f"sample: {sample}")
    alignment_file = f'/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/{sample}_primary_over_30_sorted.bam'
    alignment = list()
    alignment.append(alignment_file)
    alignment.sort()
else:
    alignment_dir = '/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/'
    alignment = list()
    for file in os.listdir(alignment_dir):
        if file.endswith("primary_over_30_sorted.bam"):
            alignment.append(os.path.join(alignment_dir, file))
    alignment.sort()



###Define the function to compute the numebr of exon-exon junctions covered by one read based on CIGAR string
def compute_num_junction_per_read(cigar_string):
    #m: number of exon-exon splice junctions
    m = cigar_string.count('N')
    return m
  
###Calculate for each bam file
max_junction = 11 #can be changed
nums = list()
df_dic = {}
# whole_cigar = ''
for fn in alignment:
    bamfile = pysam.AlignmentFile(fn, "rb") 
    for i in range(genic_gtf.shape[0]):
    # for read in bamfile.fetch('chr3', 149964904, 149970895): #PFN2, at most 2 exon-exon junctions
    # for read in bamfile.fetch('chr14', 69767112, 69772005): #SRSF5, at most 8 exon-exon junctions
        for read in bamfile.fetch(genic_gtf.iloc[i,0], genic_gtf.iloc[i,3], genic_gtf.iloc[i,4]):  #whole genome
        #     print(read.mapq)
        #     print(read)
        #     break
    #         name = read.query_name
    #         names.append(name)
    #         mq = read.mapq
    #         mqs.append(mq)
            num = compute_num_junction_per_read(read.cigarstring)
            nums.append(num)
    #         cigar = read.cigarstring
    #         whole_cigar = whole_cigar + cigar
    #     print(Counter(nums))
    tmp_counter = Counter(nums)
    # print(tmp_counter)
    tmp_num = list()
    for i in range(max_junction):   
        tmp_num.append(tmp_counter[i])
    df_dic[fn] = tmp_num    
    nums = list()

df = pd.DataFrame(df_dic).transpose()
per_df = df.div(df.sum(axis=1), axis=0)

if len(sys.argv) > 1:
    print(f"Done for {sample}")
    per_df.to_csv('/users/sparthib/retina_lrs/processed_data/exon_exon/{sample}_junction_per_read.csv')
else:
    per_df.to_csv('/users/sparthib/retina_lrs/processed_data/exon_exon/junction_per_read.csv')
    print("Done for all samples")

    
  
