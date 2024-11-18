blaze_bc_dir <- "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/03_blaze_processed/raw/high_sensitivity"

library(dplyr)
#load all files from dir ending with _putative_bc.csv
files <- list.files(blaze_bc_dir, pattern = "_whitelist.csv$", full.names = TRUE)

# Function to load a file, convert to dataframe, and add sample name column
load_and_label_file <- function(file_path) {
  # Load file into a dataframe
  df <- read.csv(file_path, header = FALSE, col.names = "barcode")
  
  # Extract the sample name from the filename (remove directory and extension)
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Add a new column called 'sample' with the sample name
  df$sample <- sample_name
  
  df$day <- ifelse(grepl("D100", sample_name), "D100",
                   ifelse(grepl("D200", sample_name), "D200", NA))
  return(df)
}

# Apply the function to each file and combine into a single dataframe
all_data <- do.call(rbind, lapply(files, load_and_label_file))

# View the combined dataframe
head(all_data)


# Filter for D100 and D200 barcodes
d100_barcodes <- all_data |>
  filter(day == "D100") |>
  pull(barcode)

d200_barcodes <- all_data |>
  filter(day == "D200") |>
  pull(barcode)

saveRDS(all_data, "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/sample_barcodes.rds")




