library("tidyverse")
library("here") 
library("rmarkdown")
library("knitr")

#formatting packages
library("kableExtra")
library("janitor")
library("scales")
library("ggpubr")
library("readr")


library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# 9-BRN3B vs EP1-BRN3B ROs

EP1_BRN3B_transcript_counts <- read_tsv("/dcs04/hicks/data/sparthib/casey/IsoQuant_output/EP1-BRN3B-RO/OUT/OUT.transcript_counts.tsv",
                               col_names=TRUE)
H9_BRN3B_transcript_counts <- read_tsv("/dcs04/hicks/data/sparthib/casey/IsoQuant_output/H9-BRN3B-RO/OUT/OUT.transcript_counts.tsv",
                            col_names=TRUE)



small_counts <- read.table("data/small_counts.txt", header = TRUE)