

#read in classification.txt file and tabulate the 6th column 
input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/sqanti3_qc"

#QC 
qc_classification <- read.table(paste0(input_dir, "/all_samples_classification.txt"),
                                header = TRUE, sep = "\t")
table(qc_classification$structural_category)


# antisense       full-splice_match                  fusion 
# 8139                  117747                     667 
# genic incomplete-splice_match              intergenic 
# 1904                    8380                    6483 
# novel_in_catalog    novel_not_in_catalog 
# 6603                   14160

# > nrow(qc_classification)
# [1] 164083

#output
filter_classification <- read.table(paste0(input_dir, "/all_samples_RulesFilter_result_classification.txt"),
                                    header = TRUE, sep = "\t")
table(filter_classification$structural_category)


#rescue 
rescue_table <- read.table(paste0(input_dir, "/post_rescue_automatic_rescue_table.tsv"),
                            header = TRUE, sep = "\t")
table(rescue_table$structural_category)

# table(rescue_table$structural_category)
# 
# full-splice_match incomplete-splice_match 
# 19093                     555 



### filter out genic-genomic, genic-intron, fusion, 


