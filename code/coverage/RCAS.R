library(RCAS)

queryRegions <- importBed("/dcs04/hicks/data/sparthib/retina_lrs/05b_beds/EP1-BRN3B-RO_sorted.bed")
gff <- importGtf(filePath = "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf")


runReport( queryFilePath = 'input.BED',
           gffFilePath = 'annotation.gtf', 
           genomeVersion = 'hg19',
           goAnalysis=FALSE,
           outDir = "/users/sparthib/retina_lrs/processed_data/coverage/RCAS")

