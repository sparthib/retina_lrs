### compare Isoquant single cell D100 and bulk RNA-seq data 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)


###get D100 DTU files 

DGE_DTU_DTE <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DGE_DTU_DTE.tsv")

###get D100 DTE files

D100_high_DTE <- DGE_DTU_DTE |> filter(DTE == TRUE) |> filter(DTE_qval < 0.05) |>
  filter( (condition_1 ==  "B_RO_D100" & DTE_log2FC >= 1) | 
                                                (condition_2 ==  "B_RO_D100" & DTE_log2FC < -1 )) |>
  select(isoform_id) |> distinct()

D200_high_DTE <- DGE_DTU_DTE |> filter(DTE == TRUE) |> filter(DTE_qval < 0.05) |>
  filter( (condition_1 ==  "A_RO_D200" & DTE_log2FC >= 1) | 
                                                (condition_2 ==  "A_RO_D200" & DTE_log2FC < -1)) |>
  select(isoform_id) |> distinct()


#single cell data dir 

sample_names <- c("10x_D100-EP1", "10x_D200-EP1-1",
                  "10x_D200-EP1-2" )

D100_sc_dir = paste0("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/", 
                     sample_names[1], "/OUT")

D100_sc_transcript_tpm = read.table(paste(D100_sc_dir, "OUT.transcript_model_grouped_tpm.tsv", sep="/"), 
                                    header=T, row.names=1, sep="\t", comment.char = "") 


#check for duplicates in D100_sc_transcript_tpm
duplicates = duplicated(colnames(D100_sc_transcript_tpm))


#average tpm values for each gene across all cells
D100_sc_transcript_average_tpm = rowMeans(D100_sc_transcript_tpm)

#remove values less than 10
D100_sc_transcript_average_tpm = D100_sc_transcript_average_tpm[D100_sc_transcript_average_tpm > 10]

#bulk data dir
bulk_dir = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT"
bulk_transcript_tpm = read.table(paste(bulk_dir, "OUT.transcript_model_grouped_tpm.tsv", sep="/"), 
                                    header=F, row.names=1, sep="\t")
colnames(bulk_transcript_tpm) <- c( "EP1-BRN3B-RO", "EP1-WT-D45", "EP1-WT-hRO2",
  "H9-BRN3B-RO", "H9-BRN3B-hRO2", "H9-CRX-D45", "H9-CRX-hRO2",
  "H9-FT1", "H9-FT2", "H9-hRGC1", "H9-hRGC2"
)

#remove version number from ENST names

bulk_transcript_tpm <- bulk_transcript_tpm[,1:7]

# Filter rows to keep only those with "ENST" in row names
D100_sc_transcript_average_tpm <- D100_sc_transcript_average_tpm[grep("ENST", names(D100_sc_transcript_average_tpm))]
bulk_transcript_tpm <- bulk_transcript_tpm[grep("ENST", rownames(bulk_transcript_tpm)), ]


rownames(bulk_transcript_tpm) <- gsub("\\..*", "", rownames(bulk_transcript_tpm))
names(D100_sc_transcript_average_tpm) <- gsub("\\..*", "", names(D100_sc_transcript_average_tpm))



# Identify common names between sc_transcript_average_tpm and bulk_transcript_tpm
common_names <- intersect(names(D100_sc_transcript_average_tpm), rownames(bulk_transcript_tpm))

# Subset sc_transcript_average_tpm and bulk_transcript_tpm to keep only common names
D100_sc_transcript_average_tpm <- D100_sc_transcript_average_tpm[common_names]
bulk_transcript_tpm <- bulk_transcript_tpm[common_names, ]

# Check the results
length(D100_sc_transcript_average_tpm)
dim(bulk_transcript_tpm)

D100_DTE_isoforms <- intersect(common_names, D100_high_DTE$isoform_id)
length(D100_DTE_isoforms)


#DTE 
D100_sc_transcript_average_tpm_DTE = D100_sc_transcript_average_tpm[D100_DTE_isoforms]
bulk_transcript_tpm_DTE = bulk_transcript_tpm[D100_DTE_isoforms,]



# Function to compute correlation for DTE isoforms
compute_dte_correlations <- function(sc_tpm, bulk_tpm, method = "spearman") {
  # Filter data for DTE isoforms

  # Initialize correlation vector
  correlations <- numeric(ncol(bulk_tpm))
  
  # Compute correlation for each bulk column
  for (i in seq_along(bulk_tpm )) {
    correlations[i] <- cor(sc_tpm, bulk_tpm[, i], method = method)
  }
  
  # Assign names to the correlation vector
  names(correlations) <- colnames(bulk_tpm)
  
  # Return the correlation vector
  return(correlations)
}

D100_spearman_correlations_dte <- compute_dte_correlations(D100_sc_transcript_average_tpm_DTE,
                                                      bulk_transcript_tpm_DTE,
                                                      method = "spearman")

# EP1-BRN3B-RO    EP1-WT-D45   EP1-WT-hRO2   H9-BRN3B-RO H9-BRN3B-hRO2 
# 0.3854081     0.4087311     0.4271472     0.3078165     0.4574634 
# H9-CRX-D45   H9-CRX-hRO2 
# 0.4118873     0.4473093 


############# D200 transcripts #############



D200_1_sc_dir = paste0("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/", 
                     sample_names[2], "/OUT")
D200_2_sc_dir = paste0("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/", 
                       sample_names[3], "/OUT")

D200_1_sc_transcript_tpm = read.table(paste(D200_1_sc_dir, "OUT.transcript_model_grouped_tpm.tsv", sep="/"), 
                                    header=T, row.names=1, sep="\t", comment.char = "") 
D200_2_sc_transcript_tpm = read.table(paste(D200_2_sc_dir, "OUT.transcript_model_grouped_tpm.tsv", sep="/"), 
                                      header=T, row.names=1, sep="\t", comment.char = "") 



#append colnames of sc_transcript_tpms with D200_1 and D200_2
colnames(D200_1_sc_transcript_tpm) <- paste0(colnames(D200_1_sc_transcript_tpm), "_D200_1")
colnames(D200_2_sc_transcript_tpm) <- paste0(colnames(D200_2_sc_transcript_tpm), "_D200_2")



#average tpm values for each gene across all cells
D200_1_sc_transcript_average_tpm = rowMeans(D200_1_sc_transcript_tpm)
D200_2_sc_transcript_average_tpm = rowMeans(D200_2_sc_transcript_tpm)

#remove values less than 10
D200_1_sc_transcript_average_tpm = D200_1_sc_transcript_average_tpm[D200_1_sc_transcript_average_tpm > 10]
D200_2_sc_transcript_average_tpm = D200_2_sc_transcript_average_tpm[D200_2_sc_transcript_average_tpm > 10]

# Filter rows to keep only those with "ENST" in row names
D200_1_sc_transcript_average_tpm <- D200_1_sc_transcript_average_tpm[grep("ENST", names(D200_1_sc_transcript_average_tpm))]
D200_2_sc_transcript_average_tpm <- D200_2_sc_transcript_average_tpm[grep("ENST", names(D200_2_sc_transcript_average_tpm))]

#merge D200_1 and D200_2 sc_transcript_average_tpm 

D200_sc_transcript_average_tpm = merge(D200_1_sc_transcript_average_tpm,
                                       D200_2_sc_transcript_average_tpm,
                                       by="row.names", all=TRUE)
#convert NAs to 0 in D200_sc_transcript_average_tpm


D200_sc_transcript_average_tpm[is.na(D200_sc_transcript_average_tpm)] <- 0

# average x and y 

D200_sc_transcript_average_tpm$average <- rowMeans(D200_sc_transcript_average_tpm[,2:3])
rownames(D200_sc_transcript_average_tpm) <- D200_sc_transcript_average_tpm$Row.names


D200_sc_transcript_average_tpm <- D200_sc_transcript_average_tpm |> select(average)


bulk_transcript_tpm <- bulk_transcript_tpm[grep("ENST", rownames(bulk_transcript_tpm)), ]


rownames(bulk_transcript_tpm) <- gsub("\\..*", "", rownames(bulk_transcript_tpm))
rownames(D200_sc_transcript_average_tpm) <- gsub("\\..*", "", rownames(D200_sc_transcript_average_tpm))


common_names <- intersect(rownames(D200_sc_transcript_average_tpm), rownames(bulk_transcript_tpm))

# Subset sc_transcript_average_tpm and bulk_transcript_tpm to keep only common names
D200_sc_transcript_average_tpm <- D200_sc_transcript_average_tpm |> filter(rownames(D200_sc_transcript_average_tpm) 
                                                                           %in% common_names)
bulk_transcript_tpm <- bulk_transcript_tpm[common_names, ]

# Check the results
nrow(D200_sc_transcript_average_tpm)
dim(bulk_transcript_tpm)

D200_DTE_isoforms <- intersect(common_names, D200_high_DTE$isoform_id)


#DTE 
D200_sc_transcript_average_tpm_DTE = D200_sc_transcript_average_tpm |> 
  filter(rownames(D200_sc_transcript_average_tpm) %in% D200_DTE_isoforms)
                                                                              
bulk_transcript_tpm_DTE = bulk_transcript_tpm[D200_DTE_isoforms,]

D200_spearman_correlations_dte <- compute_dte_correlations(D200_sc_transcript_average_tpm_DTE$average,
                                                      bulk_transcript_tpm_DTE,
                                                      method = "spearman")
##D200 vs bulk spearman correlations for DTE isoforms
# EP1-BRN3B-RO    EP1-WT-D45   EP1-WT-hRO2   H9-BRN3B-RO H9-BRN3B-hRO2 
# 0.5040668     0.1772044     0.4378143     0.5414833     0.4159681 
# H9-CRX-D45   H9-CRX-hRO2 
# 0.1642089     0.4374109 

##D100 vs bulk spearman correlations for DTE isoforms
# EP1-BRN3B-RO    EP1-WT-D45   EP1-WT-hRO2   H9-BRN3B-RO H9-BRN3B-hRO2 
# 0.3854081     0.4087311     0.4271472     0.3078165     0.4574634 
# H9-CRX-D45   H9-CRX-hRO2 
# 0.4118873     0.4473093 


# Create data frame for D200 and D100 correlations

library(pheatmap)
correlations_df <- data.frame(
  Sample = c("EP1-BRN3B-RO", "EP1-WT-D45", "EP1-WT-hRO2", "H9-BRN3B-RO", "H9-BRN3B-hRO2", "H9-CRX-D45", "H9-CRX-hRO2"),
  D200_sc = c(0.5040668, 0.1772044, 0.4378143, 0.5414833, 0.4159681, 0.1642089, 0.4374109),
  D100_sc = c(0.3854081, 0.4087311, 0.4271472, 0.3078165, 0.4574634, 0.4118873, 0.4473093)
)

ordered_samples <- c("EP1-WT-D45", "H9-CRX-D45", "EP1-WT-hRO2", "H9-CRX-hRO2",  "H9-BRN3B-hRO2",
                     "EP1-BRN3B-RO", "H9-BRN3B-RO")
correlations_df <- correlations_df %>%
  filter(Sample %in% ordered_samples) %>%
  arrange(match(Sample, ordered_samples))

# Convert to matrix format for pheatmap
correlations_matrix <- as.matrix(correlations_df[, -1])  # Remove the 'Sample' column
rownames(correlations_matrix) <- correlations_df$Sample

# Convert data to long format for ggplot
correlations_long <- pivot_longer(correlations_df, cols = starts_with("D"), 
                                  names_to = "Condition", values_to = "Spearman_Correlation")

pdf("/users/sparthib/retina_lrs/single_cell_plots/sc_bulk_tpm_corr.pdf")
p <- pheatmap(correlations_matrix, 
              cluster_rows = FALSE,  # Preserve the specified order
              cluster_cols = FALSE,  # No clustering for columns
              main = "Spearman Correlations for DTE Isoforms",
              color = colorRampPalette(c("red", "white", "blue"))(50),
              display_numbers = TRUE,  # Optional: show correlation values on the heatmap
              fontsize_row = 10,
              fontsize_col = 10)
print(p)
dev.off()







