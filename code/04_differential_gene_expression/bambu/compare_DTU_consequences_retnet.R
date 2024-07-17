library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)
library(sessioninfo)
library(rtracklayer)
library(edgeR)
library(tidyr)
library(dplyr)
library(readxl)
library(here)
# library(ggVennDiagram)
library('BSgenome.Hsapiens.UCSC.hg38')

SwitchList_part2 <- readRDS( "./processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/SwitchList_part2.rds")

head(SwitchList_part2$AlternativeSplicingAnalysis)


colnames(head(SwitchList_part2$isoformFeatures))

# pre-mature termination codons
SwitchList_part2$isoformFeatures |> dplyr::filter(abs(dIF) >= 0.1) |>
  filter(isoform_switch_q_value < 0.05) |> select(PTC) |> table()

# FALSE  TRUE 
# 1096    67 

#coding potential value on scale of 0 to 1 
#cutoff suggested for humans is 0.725

SwitchList_part2$isoformFeatures |> dplyr::filter(abs(dIF) >= 0.1) |>
  filter(isoform_switch_q_value < 0.05) |> filter(codingPotentialValue > 0.725)  |> nrow()
# 1163

SwitchList_part2$isoformFeatures |> dplyr::filter(abs(dIF) >= 0.1) |>
  filter(isoform_switch_q_value < 0.05) |> select(switchConsequencesGene) |> table()
# FALSE  TRUE 
# 447   849 

SwitchList_part2$isoformFeatures |> dplyr::filter(abs(dIF) >= 0.1) |>
  filter(isoform_switch_q_value < 0.05) |> select(domain_identified) |> table()
# domain_identified
# no  yes 
# 113 1050 

SwitchList_part2$isoformFeatures |> dplyr::filter(abs(dIF) >= 0.1) |>
  filter(isoform_switch_q_value < 0.05) |> select(signal_peptide_identified) |> table()
# no  yes 
# 1055  108 

## isoforms that have PTC = FALSE, codingPotentialValue > 0.725, 
## switchConsequencesGene = TRUE, domain_identified = yes, signal_peptide_identified = yes

retn <- SwitchList_part2$isoformFeatures |> dplyr::filter(abs(dIF) >= 0.1) |>
  filter(isoform_switch_q_value < 0.05) |> filter(PTC == FALSE) |> filter(codingPotentialValue > 0.725) |>
  filter(switchConsequencesGene == TRUE) |> filter(domain_identified == "yes") |> filter(signal_peptide_identified == "yes")|>
  select(isoform_id, ensembl_gene_name)



RetNet_gene_list <- read_excel(here("raw_data", "RetNet.xlsx"),
                               sheet = "genes_and_locations")



intersect(retn$ensembl_gene_name, RetNet_gene_list$Symbol) 
"RIMS1"  "COL9A1" "COL2A1"

genes_and_diseases <- read_excel(here("raw_data", "RetNet.xlsx"),
                                 sheet = "diseases_and_genes")
colnames(genes_and_diseases) <- c("disease_category", "mapped_loci",
                                  "mapped_and_identified_genes")

results_list <- list()

# Loop through each gene
for (gene in intersect(retn$ensembl_gene_name, RetNet_gene_list$Symbol)) {
  # Filter the data
  result <- genes_and_diseases |>
    dplyr::filter(str_detect(mapped_and_identified_genes, gene))
  
  # Append results to the list
  results_list[[gene]] <- result$disease_category
}

# Convert the list to a dataframe
results_df <- data.frame(
  Gene = names(results_list),
  Disease_Category = I(results_list),
  stringsAsFactors = FALSE
)

# If you want to view the dataframe, you can use:
print(results_df)

# [1] "Cone or cone-rod dystrophy, autosomal dominant"
# 
# [[2]]
# [1] "Syndromic/systemic diseases with retinopathy, autosomal recessive"
# 
# [[3]]
# [1] "Syndromic/systemic diseases with retinopathy, autosomal dominant"
