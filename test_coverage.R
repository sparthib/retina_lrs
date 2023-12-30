# remotes::install_github("showteeth/ggcoverage")


library("rtracklayer")
library("ggcoverage")
library("ggpattern")
library("sessioninfo")

# test_bam_path <- "/dcs04/hicks/data/sparthib/casey/bams_3/H9-BRN3B-RO/H9-BRN3B-RO_sorted.bam"
sample.meta <- data.frame(
  SampleName = "DG-WT-hRGC_sorted", 
  Type = "RGC_rep1",       
  Group = "RGC")


track.df = LoadTrackFile(track.file = "/dcs04/hicks/data/sparthib/casey/bams/DG-WT-hRGC/DG-WT-hRGC_sorted.bam",
                         format = "bam", norm.method = "None",
                         meta.info = sample.meta,
                         region="chr21:1-46,709,983")

track.df

gtf.file = "/dcs04/hicks/data/sparthib/GENCODE_GTF.gtf"
gtf.gr = rtracklayer::import.gff(con = gtf.file, format = 'gtf')




mark.region=data.frame(start=5022069,
                       end=5022522,
                       label="M1")
# check data
mark.region


basic.coverage = ggcoverage(data = track.df, color = "auto", range.position = "out",
                            plot.type = "facet", 
                            mark.region = NULL)

pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/basic.coverage.pdf")

print(basic.coverage)
dev.off()



##########

# load metadata
meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
sample.meta = read.csv(meta.file)
sample.meta
#>        SampleName    Type Group
#> 1 ERR127302_chr14 KO_rep1    KO
#> 2 ERR127303_chr14 KO_rep2    KO
#> 3 ERR127306_chr14 WT_rep1    WT
#> 4 ERR127307_chr14 WT_rep2    WT



track.folder = system.file("extdata", "RNA-seq", package = "ggcoverage")
# load bigwig file
track.df = LoadTrackFile(track.folder = track.folder, format = "bw",
                         region = "chr14:21,677,306-21,737,601", extend = 2000,
                         meta.info = sample.meta)
# check data
head(track.df)

mark.region=data.frame(start=c(21678900,21732001,21737590),
                       end=c(21679900,21732400,21737650),
                       label=c("M1", "M2", "M3"))
# check data
mark.region


gtf.file = system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
gtf.gr = rtracklayer::import.gff(con = gtf.file, format = 'gtf')


pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/basic.coverage.pdf")
basic.coverage = ggcoverage(data = track.df, color = "auto", 
                            plot.type = "joint", facet.key = "Type", group.key = "Type",
                            mark.region = mark.region, range.position = "out")
print(basic.coverage)
dev.off()




##### LOCAL ######


