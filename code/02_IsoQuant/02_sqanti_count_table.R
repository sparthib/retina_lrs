library("readr")
library("here")
library("dplyr")
library("sessioninfo")

here()
### LOAD READ BARCODE DATA ###

data_paths <- read_csv(file = here("raw_data", "data_paths.csv") , 
                     col_names = T)

sample_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

class(data_paths$Index)

read_assignment_path <- paste0("/dcs04/hicks/data/sparthib/casey/IsoQuant_output/",
                               data_paths$sample_name[sample_num],
                               "/OUT/OUT.read_assignments.tsv")

## LOAD ISOQUANT OUTPUT ###
isoquant_read_assignments <- readr::read_tsv(read_assignment_path,
                                             skip = 2,
                                   col_names = TRUE)

print("assignment type")
print(table(isoquant_read_assignments$assignment_type))
print(table(isoquant_read_assignments$chrom))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


