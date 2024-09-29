# Assuming you have 4 matrices: matrix1, matrix2, matrix3, matrix4
library(readr)
library(tibble)
#Group Isoquant_matrices 
sample_names <- c( "10x_D100-EP1_A1",
                   "10x_D100-EP1_A2",
                   "10x_D100-EP1_B1",
                   "10x_D100-EP1_B2",
                   "10x_D200-EP1-1_A1",
                   "10x_D200-EP1-1_A2",
                   "10x_D200-EP1-1_B1",
                   "10x_D200-EP1-1_B2",
                   "10x_D200-EP1-2_A1",
                   "10x_D200-EP1-2_A2",
                   "10x_D200-EP1-2_B1",
                   "10x_D200-EP1-2_B2")


#read a csv.gz file

load_matrix <- function(sample) {
  # Replace invalid characters (e.g., '-') with '_'
  sanitized_sample <- gsub("-", "_", sample)
  matrix_file <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant", 
                           paste0(sample, "/OUT/OUT.transcript_model_grouped_tpm.tsv"))
  
  # Check if the file exists before attempting to read it
  if (!file.exists(matrix_file)) {
    stop("File not found: ", matrix_file)
  }
  
  # Read the gzipped CSV file as a data frame
  df <- read_tsv(matrix_file, col_names = TRUE)
  
  df <- column_to_rownames(df, var = colnames(df)[1])
  # Convert all values in the data frame to numeric
  df[] <- lapply(df, as.numeric)  # Convert each column to numeric
  
  # Convert the data frame back to a matrix
  mat <- as.matrix(df)
  
  # Assign the matrix to a variable in the global environment
  assign(paste0("matrix_", sanitized_sample), mat, envir = .GlobalEnv)
}

# Example usage
# load_matrix("sample_name")


# Example usage
# load_matrix("sample_name")


# Example usage
# load_matrix("sample_name")



mat_list <- sample_names |> lapply(load_matrix)

rownames(mat_list[[1]])

mat_list[[1]][1:5, 1:5]

output_dir <- "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/transcriptome/isoquant/merged_matrices"

# Function to merge and sum counts matrices
merge_counts_matrices <- function(matrices, fill_value = 0, output_name) {
  # Capture all input matrices into a list
  
  # Get the union of row and column names across all matrices
  all_rows <- Reduce(union, lapply(matrices, rownames))
  all_cols <- Reduce(union, lapply(matrices, colnames))
  
  # Initialize a list to store filled matrices
  filled_matrices <- lapply(matrices, function(mat) {
    # Create an empty matrix with all possible row and column names, filled with zeros (or specified fill_value)
    filled_mat <- matrix(fill_value, nrow=length(all_rows), ncol=length(all_cols), dimnames=list(all_rows, all_cols))
    
    # Copy values from the original matrix to the new one
    filled_mat[rownames(mat), colnames(mat)] <- mat
    return(filled_mat)
  })
  
  # Sum all the filled matrices element-wise
  merged_matrix <- Reduce(`+`, filled_matrices)
  
  print("finished merging matrices")
  saveRDS( merged_matrix, file.path(output_dir, paste0(output_name, ".rds")))
  # Return the merged matrix
  # return(merged_matrix)
}

# Example usage with matrix1, matrix2, matrix3, matrix4
merge_counts_matrices(mat_list[1:4], fill_value = 0, output_name = "mat_10x_D100_EP1")
merge_counts_matrices(mat_list[5:8], fill_value = 0, output_name = "mat_10x_D200_EP1_1")
merge_counts_matrices(mat_list[9:12], fill_value = 0, output_name = "mat_10x_D200_EP1_2")


mat_10x_D100_EP1 <- readRDS(file.path(output_dir, "mat_10x_D100_EP1.rds"))
mat_10x_D200_EP1_1 <- readRDS(file.path(output_dir, "mat_10x_D200_EP1_1.rds"))
mat_10x_D200_EP1_2 <- readRDS(file.path(output_dir, "mat_10x_D200_EP1_2.rds"))

dim(mat_10x_D100_EP1)
dim(mat_10x_D200_EP1_1)
dim(mat_10x_D200_EP1_2)






## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



