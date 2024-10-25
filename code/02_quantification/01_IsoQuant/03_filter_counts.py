import pandas as pd

# https://github.com/RCHENLAB/LR_scRNA-seq_manuscript/blob/main/LR_analysis/4isoform_filtering.py




# Read in the cell class by isoform count matrix and SQANTIs isoform classification
counts_file="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT/OUT.transcript_model_grouped_counts.tsv"

classification_text="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/sqanti3_qc/all_samples_MLresult_classification.txt"

rescued_counts_file="/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT/final_rescued_counts.tsv"
#read tsv on pandas 
list1 = pd.read_csv(counts_file, sep = "\t", encoding = "ISO-8859-1", header = (0))
list2 = pd.read_csv(classification_text, sep = "\t", encoding = "ISO-8859-1", header = (0))
list1.columns.values[0] = 'id'
list2.columns.values[0] = 'id'

df = list1.merge(list2, how='inner', on=["id"])
# numbers of rows of df
print(df.shape[0])
#85584

#get column names of df
cols = df.columns.tolist()
# Filter out isoforms classified as artifacts by SQANTI3
df_filtered = df.loc[(df['filter_result'] == 'Isoform')]
df_artifacts = df.loc[(df['filter_result'] == 'Artifact')]
# numbers of rows of df_filtered
print(df_filtered.shape[0])
# 52380
# Write the output to a file
df.to_csv(rescued_counts_file, index = False, sep = "\t")

#get structural categories of isoforms from SQANTI3 classification
df_filtered['structural_category'].value_counts()
# structural_category
# full-splice_match          31958
# intergenic                  5195
# antisense                   5194
# incomplete-splice_match     4034
# novel_not_in_catalog        2687
# novel_in_catalog            2233
# genic                        786
# fusion                       293

df_artifacts['structural_category'].value_counts()
# >>> df_artifacts['structural_category'].value_counts()
# structural_category
# novel_not_in_catalog       11467
# full-splice_match           7302
# novel_in_catalog            4366
# incomplete-splice_match     4344
# antisense                   2945
# intergenic                  1288
# genic                       1118
# fusion                       374



