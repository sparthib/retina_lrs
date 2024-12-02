library("rtracklayer")

sqanti_gtf = rtracklayer::import("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/sqanti3_qc/sqanti3_qc/FT_vs_RGC_corrected.gtf.cds.gff")
tx.sqanti = sqanti_gtf |> as_tibble() |> filter(type == "transcript")

