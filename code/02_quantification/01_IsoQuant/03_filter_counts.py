import pandas as pd

# Read in the cell class by isoform count matrix and SQANTIs isoform classification
counts_file="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT/OUT.transcript_model_grouped_counts.tsv"
classification_text="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/sqanti3_qc/all_samples_classification.txt"
rescued_counts_file="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT/final_rescued_counts.tsv"
#read tsv on pandas 
list1 = pd.read_csv(counts_file, sep = "\t", encoding = "ISO-8859-1", header = (0))
list2 = pd.read_csv(classification_text, sep = "\t", encoding = "ISO-8859-1", header = (0))
list1.columns.values[0] = 'id'
list2.columns.values[0] = 'id'

df = list1.merge(list2, how='inner', on=["id"])

# Write the output to a file
df.to_csv(rescued_counts_file, index = False, sep = "\t")
