library("readr")
library("here")
library("dplyr")
library("sessioninfo")
library("tidyr")

here()
### LOAD READ BARCODE DATA ###

data_paths <- read_csv(file = here("raw_data", "data_paths.csv") , 
                     col_names = T)

sample_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

print(data_paths$sample_name[sample_num])

read_assignment_path <- paste0("/dcs04/hicks/data/sparthib/casey/IsoQuant_output/",
                               data_paths$sample_name[sample_num],
                               "/OUT/OUT.read_assignments.tsv")

## LOAD ISOQUANT OUTPUT ###
isoquant_read_assignments <- readr::read_tsv(read_assignment_path,
                                   col_names = TRUE)

print("assignment type")
print(table(isoquant_read_assignments$chr))

add_info <- separate_wider_delim(isoquant_read_assignments, 
                           cols = additional_info, delim = ";", 
                           names = c("polyA", "Canonical", "Classification"))

print(table(add_info$Classification))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


