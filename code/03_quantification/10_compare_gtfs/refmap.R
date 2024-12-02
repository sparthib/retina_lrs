######

refmap_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT/"
refmap <- "bambu_ref_isoquant_ext.OUT.extended_annotation.gtf.refmap"

refmap_path <- paste(refmap_dir, refmap, sep="")

refmap <- read.table(refmap_path, header=TRUE, sep="\t")

refmap[1:5,1:4]

#filter out ref_id where isoform doesn't start with ENST

bambu_novel_refmap <- refmap[!grepl("^ENST", refmap$ref_id), ]
nrow(bambu_novel_refmap)

bambu_novel_refmap <- refmap[!grepl("^ENST", refmap$ref_id) & refmap$class_code == "=", ]

# nrow(bambu_novel_refmap)
#599 

#write out the refmap
# /dcs04/hicks/data/sparthib/retina_lrs/06_quantification/gffcompare
output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/"
write.table(bambu_novel_refmap, file=paste(output_dir, "bambu_isoform_novel_refmap.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

