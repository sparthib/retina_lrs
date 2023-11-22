library(optparse)
library(ShortRead)

parser <- OptionParser()
parser <- add_option(parser, 
                     opt_str = c("-i", "--input"), 
                     type = "character",
                     dest = 'input.file',
                     help="Input file or directory (required). fastq file")

parser <- add_option(parser, 
                     opt_str = c("-o", "--output"), 
                     type = "character",
                     dest = 'output.file',
                     help="output file path")



opt = parse_args(parser)
if (length(opt$input.file)==1) {
  input.file = opt$input.file
} else {
  stop("Input file parameter must be supplied via -i or --input")
}



f <- FastqStreamer(input.file,readerBlockSize=1000) 

# we set up a while loop to call yield() function to
# go through the file
while(length(fq <- yield(f))) {
  
  # remove reads where all quality scores are < 20 
  # get quality scores per base as a matrix
  qPerBase = as(quality(fq), "matrix")
  
  # get number of bases per read that have Q score < 20
  qcount = rowSums( qPerBase <= 20) 
  
  # write fastq file with mode="a", so every new block
  # is written out to the same file
  writeFastq(fq[qcount == 0], 
             paste(fastqFile, "Qfiltered", sep="_"), 
             mode="a")
}