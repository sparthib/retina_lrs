import pysam
from collections import Counter
import pandas as pd
import numpy as np
import os
import pyranges as pr
import sys

#from https://github.molgen.mpg.de/MayerGroup/NGN3_paper_code/blob/main/ONT_seq_analysis/long_short_comparison/junction_per_read.py
# Load the GTF file
# Load the GTF file
gtf_fn = "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
gtf = pd.read_csv(gtf_fn, sep='\t', header=None, comment='#', 
                  names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])


curr_gene = ''
goi = list()
curr_num_exon = 0
for i in range(gtf.shape[0]):
    if gtf.iloc[i, 2] == 'gene' and 'protein_coding' in gtf.iloc[i, 8]:  # Adjust to check for 'protein_coding'
        if curr_num_exon > 1:
            goi.append(curr_gene)
        curr_gene = gtf.iloc[i, 8].split(';')[0]  # Get the gene ID
        curr_num_exon = 0
    if gtf.iloc[i, 2] == 'exon' and curr_gene in gtf.iloc[i, 8]:
        curr_num_exon += 1

genic_gtf = gtf[gtf['feature'] == 'gene']
genic_gtf = genic_gtf[genic_gtf['attribute'].str.contains('|'.join(goi))]
genic_gtf = genic_gtf.iloc[:, 0:14]

###Get bam file from command line arg
if len(sys.argv) > 1:
    sample = sys.argv[1]
    print(f"sample: {sample}")

    alignment_file = f'/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/{sample}_primary_over_30_chr_only_sorted.bam'
    alignment = [alignment_file]
    alignment.sort()
else:
    alignment_dir = '/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality/'
    alignment = []
    for file in os.listdir(alignment_dir):
        if file.endswith("primary_over_30_chr_only_sorted.bam"):
            alignment.append(os.path.join(alignment_dir, file))
    alignment.sort()



###Define the function to compute the numebr of exon-exon junctions covered by one read based on CIGAR string
def compute_num_junction_per_read(cigar_string):
    return cigar_string.count('N')

# Calculate for each bam file
max_junction = 11  # Can be changed
nums = []
df_dic = {}

for fn in alignment:
    bamfile = pysam.AlignmentFile(fn, "rb")
    for i in range(genic_gtf.shape[0]):
        contig = genic_gtf.iloc[i, 0]
        if contig in bamfile.references:  # Check if the contig exists in the BAM file
            for read in bamfile.fetch(contig, genic_gtf.iloc[i, 3], genic_gtf.iloc[i, 4]):
                num = compute_num_junction_per_read(read.cigarstring)
                nums.append(num)
        
    tmp_counter = Counter(nums)
    tmp_num = [tmp_counter.get(i, 0) for i in range(max_junction)]  # Handle missing keys
    df_dic[fn] = tmp_num
    nums = []

df = pd.DataFrame(df_dic).transpose()
per_df = df.div(df.sum(axis=1), axis=0)

# Check if a sample was passed as a command line argument
if len(sys.argv) > 1:
    sample = sys.argv[1]
    print(f"Done for {sample}")
    output_file = f'/users/sparthib/retina_lrs/processed_data/exon_exon/{sample}_junction_per_read.csv'
else:
    print("Done for all samples")
    output_file = '/users/sparthib/retina_lrs/processed_data/exon_exon/junction_per_read.csv'

# Save the DataFrame to the appropriate file
per_df.to_csv(output_file)


    
  
