library('RCAS')
library('BSgenome.Hsapiens.UCSC.hg38')
# queryRegions <- importBed("/dcs04/hicks/data/sparthib/retina_lrs/05b_beds/EP1-BRN3B-RO_sorted.bed")
# gff <- importGtf(filePath = "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf")


runReport( queryFilePath = "/dcs04/hicks/data/sparthib/retina_lrs/05b_beds/EP1-BRN3B-RO_sorted.bed",
           gffFilePath = "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf", 
           genomeVersion = 'hg38',
           goAnalysis=FALSE,
           outDir = "/users/sparthib/retina_lrs/processed_data/coverage/RCAS")

