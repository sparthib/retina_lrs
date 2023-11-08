library("readr")
library("here")
library("dplyr")
library("sessioninfo")

here()
### LOAD READ BARCODE DATA ###

data_paths <- read_csv(file = here("raw_data", "data_paths.csv") , 
                     col_names = T)
View(data_paths)

sample_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

class(data_paths$Index)


## LOAD ISOQUANT OUTPUT ###
isoquant_read_assignments <- readr::read_tsv(here("dcs04",
                                                  "hicks",
                                                  "data",
                                                  "sparthib",
                                                  "casey",
                                                  "IsoQuant_output",
                                                  data_paths$sample_name[sample_num],
                                                  "OUT",
                                                  "OUT.read_assignments.tsv"),
                                             skip = 2,
                                   col_names = TRUE)

print("assignment type")
print(table(isoquant_read_assignments$assignment_type))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


