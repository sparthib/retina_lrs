library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/GAlignments"


#if file doesn't exit,do the following steps 
if (!file.exists(file.path(output_dir, "all_samples_transcript_info.rds"))) {
  #list all .rds files in the output directory
  rds_files <- list.files(output_dir, pattern = "_transcript_info.rds$", full.names = TRUE)
  
  #read all files in rds_files
  transcript_info <- lapply(rds_files, readRDS)
  
  #add sample name to each element of the list
  sample.names <- gsub("_transcript_info.rds", "", basename(rds_files))
  sample_runs <- c("run_2", "run_3", "run_1", "run_3", 
                   "run_4", "run_4", "run_4", "run_4", 
                   "run_1", "run_3", "run_3")
  #append sample run to sample name
  sample.names <- paste(sample.names, sample_runs, sep = "_")
  transcript_info <- Map(function(df, sample) {
    df$sample_name <- sample
    return(df)
  }, transcript_info, sample.names)
  
  #add run number to each element of the list
  transcript_info <- Map(function(df, run) {
    df$run_number <- run
    return(df)
  }, transcript_info, sample_runs)
  
  #combine all dataframes into one
  transcript_info <- do.call(rbind, transcript_info)


  #save transcript_info as rds file
  saveRDS(transcript_info, file = file.path(output_dir, "all_samples_transcript_info.rds"))
} else {
  #read transcript_info from rds file
  transcript_info <- readRDS(file = file.path(output_dir, "all_samples_transcript_info.rds")) } 


# 230912_Casey	hRGCs
# 230912_Casey	EP1-BRN3B-ROs
# 230912_Casey	H9-BRN3B-ROs
# 230920_Casey	DG-WT-hRGCs
# 193679_Zack	193671_1; H9-CRX_ROs_D45
# 193679_Zack	193671_2; EP1-WT_ROs_D45
# 193679_Zack	193671_3; YZ-3T_hRGCs
# 193679_Zack	193671_4; YZ-15T_hRGCs
# CK_01052024	H9-hRGC_1
# CK_01052024	H9-FT_1
# CK_01052024	H9-hRGC_2
# CK_01052024	H9-FT_2


#group by GC content and average columns that have "coverage" in their name
gc_wise_average <- transcript_info |>
  group_by(gc_content_bin) |> 
  summarise(across(contains("coverage"), mean, na.rm = TRUE))

run_wise_average <- transcript_info |>
  group_by(run_number) |> 
  summarise(across(contains("coverage"), mean, na.rm = TRUE))
#group by transcript biotype and average columns that have "coverage" in their name
biotype_wise_average <- transcript_info |> 
  group_by(transcript_biotype) |> 
  summarise(across(contains("coverage"), mean, na.rm = TRUE))

length_bin_wise_average <- transcript_info |> 
  group_by(length_bin) |> 
  summarise(across(contains("coverage"), mean, na.rm = TRUE))


saveRDS(gc_wise_average, file = file.path(output_dir, "gc_wise_average.rds"))
saveRDS(run_wise_average, file = file.path(output_dir, "run_wise_average.rds"))
saveRDS(biotype_wise_average, file = file.path(output_dir, "biotype_wise_average.rds"))
saveRDS(length_bin_wise_average, file = file.path(output_dir, "length_bin_wise_average.rds"))

gc_wise_average <- readRDS(file = file.path(output_dir, "gc_wise_average.rds"))
run_wise_average <- readRDS(file = file.path(output_dir, "run_wise_average.rds"))
biotype_wise_average <- readRDS(file = file.path(output_dir, "biotype_wise_average.rds"))
length_bin_wise_average <- readRDS(file = file.path(output_dir, "length_bin_wise_average.rds"))

#sample specific for above categories 
sample_specific_gc_wise_average <- transcript_info |>
  group_by(sample_name, gc_content_bin) |> 
  summarise(across(contains("coverage"), mean, na.rm = TRUE))

sample_specific_run_wise_average <- transcript_info |>
  group_by(sample_name, run_number) |> 
  summarise(across(contains("coverage"), mean, na.rm = TRUE))

sample_specific_biotype_wise_average <- transcript_info |>
  group_by(sample_name, transcript_biotype) |> 
  summarise(across(contains("coverage"), mean, na.rm = TRUE))

sample_specific_length_bin_wise_average <- transcript_info |>
  group_by(sample_name, length_bin) |> 
  summarise(across(contains("coverage"), mean, na.rm = TRUE))

saveRDS(sample_specific_gc_wise_average, file = file.path(output_dir, "sample_specific_gc_wise_average.rds"))
saveRDS(sample_specific_biotype_wise_average, file = file.path(output_dir, "sample_specific_biotype_wise_average.rds"))
saveRDS(sample_specific_length_bin_wise_average, file = file.path(output_dir, "sample_specific_length_bin_wise_average.rds"))
saveRDS(sample_specific_run_wise_average, file = file.path(output_dir, "sample_specific_run_wise_average.rds"))

sample_specific_gc_wise_average <- readRDS(file = file.path(output_dir, "sample_specific_gc_wise_average.rds")) 
sample_specific_biotype_wise_average <- readRDS(file = file.path(output_dir, "sample_specific_biotype_wise_average.rds"))
sample_specific_length_bin_wise_average <- readRDS(file = file.path(output_dir, "sample_specific_length_bin_wise_average.rds"))
sample_specific_run_wise_average <- readRDS(file = file.path(output_dir, "sample_specific_run_wise_average.rds"))

#plot sample_specific_gc_wise_average

########
# Define the function
plot_coverage_trend <- function(data, group_col, output_dir, name, coverage_prefix = "coverage_", x_axis_values = 1:100) {
  
  # Convert group_col to symbol for dynamic evaluation
  group_col <- sym(group_col)
  
  # Check the actual column names
  coverage_columns <- grep(paste0("^", coverage_prefix), colnames(data), value = TRUE)
  
  # Ensure we have some coverage columns
  if (length(coverage_columns) == 0) {
    stop("No columns matching the coverage prefix were found in the data.")
  }
  
  # Reshape the data into long format
  df_long <- data |>
    pivot_longer(cols = all_of(coverage_columns), names_to = "coverage", values_to = "value") %>%
    group_by(!!group_col)
  
  # Create x-axis values for plotting
  df_long <- df_long |>
    mutate(x_axis = as.numeric(sub("coverage_", "", coverage)))
  
  # Initialize ggplot
  p <- ggplot(df_long, aes(x = x_axis, y = value, color = !!group_col, group = !!group_col)) +
    geom_line() +
    geom_smooth(se = FALSE, method = "loess") +
    labs(x = "5' to 3'", y = "Coverage", title = paste("Coverage by", as.character(group_col))) +
    theme_minimal()
  
  # Check if sample_name column is present and facet by sample if it is
  if ("sample_name" %in% colnames(data)) {
    p <- p + facet_wrap(~sample_name)
  }
  
  # Move legend to bottom and make it smaller for transcript biotype
  if (group_col == "transcript_biotype") {
    p <- p + theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )
  }
  
  # Create the file path
  file_path <- file.path(output_dir, paste0(name, ".pdf"))
  
  # Save the plot as PDF
  ggsave(file_path, plot = p, device = "pdf", width = 8, height = 6)
}


#######
coverage_prefix = "coverage_"

data = sample_specific_gc_wise_average
coverage_columns <- grep(paste0("^", coverage_prefix), colnames(data), value = TRUE)
df_long <- data |>
  pivot_longer(cols = all_of(coverage_columns), names_to = "coverage", values_to = "value") |>
  group_by(gc_content_bin)

#create a column for x-axis values
#if coverage = coverag_1 then x_axis = 1, if coverage = coverage_2 then x_axis = 2 and so on
df_long <- df_long |>
  mutate(x_axis = as.numeric(sub("coverage_", "", coverage)))

p <- ggplot(df_long, aes(x = x_axis, y = value, color = gc_content_bin, group = gc_content_bin)) +
  geom_line() +
  geom_smooth(se = FALSE, method = "loess") +
  labs(x = "X-axis (1 to 100)", y = "Coverage", title = paste("Coverage by GC content bin")) +
  theme_minimal() + facet_wrap(~sample_name)
file_path <- file.path(output_plot_dir, paste0("sample_specific_gc_wise_average", ".pdf"))
ggsave(file_path, plot = p, device = "pdf", width = 8, height = 6)


# Define the function






output_plot_dir <- "/users/sparthib/retina_lrs/plots/coverage/flames_cov_plots/"

plot_coverage_trend(run_wise_average, group_col = "run_number", output_plot_dir, "run_wise_average")
plot_coverage_trend(gc_wise_average, group_col = "gc_content_bin",  output_plot_dir,  "gc_wise_average")
plot_coverage_trend(biotype_wise_average,group_col = "transcript_biotype", output_plot_dir, "biotype_wise_average")
plot_coverage_trend(length_bin_wise_average, group_col = "length_bin", output_plot_dir, "length_bin_wise_average")
plot_coverage_trend(sample_specific_gc_wise_average, group_col = "gc_content_bin", output_plot_dir, "sample_specific_gc_wise_average")
plot_coverage_trend(sample_specific_biotype_wise_average, group_col = "transcript_biotype", output_plot_dir, "sample_specific_biotype_wise_average")
plot_coverage_trend(sample_specific_length_bin_wise_average, group_col = "length_bin", output_plot_dir, "sample_specific_length_bin_wise_average")
plot_coverage_trend(sample_specific_run_wise_average, group_col = "run_number", output_plot_dir, "sample_specific_run_wise_average")





