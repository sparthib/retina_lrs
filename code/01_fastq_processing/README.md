`01_concat_fastqs` concats fastqs output from all flowcells into one fastq per sample.

`02_fastq_qc`      filters reads based on mean read Qscore (threshold =7) and 
                   read length (minimum = 50bp) using nanofilt

`03_restrander`    explored tool for orienting reads

`04_minIONQC`      minIONQC is a tool for getting read summary stats and plots from guppy summary 
                  text file
`05_fastq_to_bam` script for aligning reads using minimap2
