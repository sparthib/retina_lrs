library(ggplot2)
library(stringr)
library(sessioninfo)

sample <- commandArgs(trailingOnly = TRUE)[1]
dir.create(paste0("/users/sparthib/retina_lrs/plots/longshot/", sample), showWarnings = FALSE)
output_dir <- paste0("/users/sparthib/retina_lrs/plots/longshot/", sample)

#Read Length distribution 
seq_dir <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/01_longshot_vcfs/",sample, "/")
HP1_read_length_distribution <- paste0(seq_dir, sample, ".HP1.read_length_distribution.txt")
HP2_read_length_distribution <- paste0(seq_dir, sample, ".HP2.read_length_distribution.txt")
noHP_read_length_distribution <- paste0(seq_dir, sample, ".noHP.read_length_distribution.txt")

HP1_read_length_distribution <- read.table(HP1_read_length_distribution, header = FALSE)
HP2_read_length_distribution <- read.table(HP2_read_length_distribution, header = FALSE)
noHP_read_length_distribution <- read.table(noHP_read_length_distribution, header = FALSE)

colnames(HP1_read_length_distribution) <- c("read_length", "count")
colnames(HP2_read_length_distribution) <- c("read_length", "count")
colnames(noHP_read_length_distribution) <- c("read_length", "count")

HP1_read_length_distribution$HP <- "HP1"
HP2_read_length_distribution$HP <- "HP2"
noHP_read_length_distribution$HP <- "noHP"

read_length_distribution <- rbind(HP1_read_length_distribution, HP2_read_length_distribution, noHP_read_length_distribution)

# Plot histogram of read length distribution using ggplot2 based on HP
pdf(paste0(output_dir, "/", sample, "_read_length_distribution.pdf"))
p <- ggplot(read_length_distribution, aes(x = read_length, y = count, fill = HP)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Read Length Distribution", x = "Read Length", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "top")
print(p)
dev.off()

print("plotted read length distribution")

#chr distribution
HP1_chr_distribution <- paste0(seq_dir, sample, ".HP1.read_chr_distribution.txt")
HP2_chr_distribution <- paste0(seq_dir, sample, ".HP2.read_chr_distribution.txt")
noHP_chr_distribution <- paste0(seq_dir, sample, ".noHP.read_chr_distribution.txt")

HP1_chr_distribution <- read.table(HP1_chr_distribution, header = FALSE)
HP2_chr_distribution <- read.table(HP2_chr_distribution, header = FALSE)
noHP_chr_distribution <- read.table(noHP_chr_distribution, header = FALSE)

colnames(HP1_chr_distribution) <- c("chr", "count")
colnames(HP2_chr_distribution) <- c("chr", "count")
colnames(noHP_chr_distribution) <- c("chr", "count")

HP1_chr_distribution$HP <- "HP1"
HP2_chr_distribution$HP <- "HP2"
noHP_chr_distribution$HP <- "noHP"

chr_distribution <- rbind(HP1_chr_distribution, HP2_chr_distribution, noHP_chr_distribution)

# Plot barplot of chr distribution using ggplot2 based on HP
pdf(paste0(output_dir, "/", sample, "_chr_distribution.pdf"))
p <- ggplot(chr_distribution, aes(x = chr, y = count, fill = HP)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Chromosome Distribution", x = "Chromosome", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "top")
print(p)
dev.off()

print("plotted chr distribution")

#GC content plot 
HP1_seq <- paste0(seq_dir, sample, ".HP1.sequences.txt")
HP2_seq <- paste0(seq_dir, sample, ".HP2.sequences.txt")
noHP_seq <- paste0(seq_dir, sample, ".noHP.sequences.txt")

HP1_sequences <- readLines(HP1_seq)
HP2_sequences <- readLines(HP2_seq)
noHP_sequences <- readLines(noHP_seq)

gc_content <- function(seq) {
  gc <- sum(str_count(seq, "G")) + sum(str_count(seq, "C"))
  return( gc * 100 / nchar(seq) )
}

HP1_gc_values <- sapply(HP1_sequences, gc_content)
names(HP1_gc_values) <- "HP1"
HP2_gc_values <- sapply(HP2_sequences, gc_content)
names(HP2_gc_values) <- "HP2"
noHP_gc_values <- sapply(noHP_sequences, gc_content)
names(noHP_gc_values) <- "noHP"

HP1_gc_df <- data.frame(GC_Content = HP1_gc_values, HP = "HP1")
HP2_gc_df <- data.frame(GC_Content = HP2_gc_values, HP = "HP2")
noHP_gc_df <- data.frame(GC_Content = noHP_gc_values, HP = "noHP")

gc_df <- rbind(HP1_gc_df, HP2_gc_df, noHP_gc_df)

# Plot histogram of GC content distribution using ggplot2 based on HP
pdf(paste0(output_dir, "/", sample, "_GC_content_distribution.pdf"))
p <- ggplot(gc_df, aes(x = GC_Content, fill = HP)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  labs(title = "GC Content Distribution", x = "GC Content", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "top")
print(p)
dev.off()

print("plotted GC content distribution")

sessionInfo()
