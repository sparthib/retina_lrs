library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(IsoformSwitchAnalyzeR)


SwitchList_part2 <- readRDS("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/rds/SwitchList_part2.rds")


# SwitchList_part2$isoformFeatures              SwitchList_part2$isoformRepExpression         SwitchList_part2$aaSequence
# SwitchList_part2$exons                        SwitchList_part2$runInfo                      SwitchList_part2$AlternativeSplicingAnalysis
# SwitchList_part2$conditions                   SwitchList_part2$orfAnalysis                  SwitchList_part2$domainAnalysis
# SwitchList_part2$designMatrix                 SwitchList_part2$isoformRepIF                 SwitchList_part2$signalPeptideAnalysis
# SwitchList_part2$sourceId                              
# SwitchList_part2$isoformCountMatrix                 


SwitchList_part2$switchConsequence |> dplyr::select(featureCompared)  |> 
  head()

# 1   intron_retention
# 2   coding_potential
# 3 ORF_seq_similarity
# 4         NMD_status
# 5 domains_identified
# 6     domain_isotype

saveRDS()
mee <- SwitchList_part2$AlternativeSplicingAnalysis |>  filter(MEE > 0 ) 

mee <- merge(mee,SwitchList_part2$isoformSwitchAnalysis, by = "isoform_id")
write_tsv(mee, "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/mee.tsv")




