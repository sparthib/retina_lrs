output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/GAlignments/"
#get all files that are transcript_info.rds in output dir
files <- list.files(output_dir, pattern = "_transcript_info.rds", full.names = TRUE)

#read all files into a list
transcript_info_list <- lapply(files, readRDS)

#for each transcript_info, add sample name
for (i in 1:length(transcript_info_list)) {
  transcript_info_list[[i]]$sample <- basename(files[i]) |> str_remove("_transcript_info.rds")
}

#combine all transcript_info into one dataframe
transcript_info <- do.call(rbind, transcript_info_list)

#group by transcript and average all columns that have "coverage" in the name
transcript_info <- transcript_info |> 
  dplyr::group_by(isoform) |> 
  dplyr::summarise(across(contains("coverage"), mean, na.rm = TRUE)) |> 
  dplyr::ungroup()

