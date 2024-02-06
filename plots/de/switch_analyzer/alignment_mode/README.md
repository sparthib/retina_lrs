This folder Contains plots run by IsoformSwitchAnalyzer on input 
from salmon run in alignment mode with `--ont` flag on. 
Accompanying `tsv` files in `processed_data/dtu/IsoformSwitchAnalyzer/salmon_alignment_mode`

`switch_vs_degs.pdf`: Pairwise Switch vs DEG plots, genes that have |log2fold change| < 2 and |dIF| > 0.5 are annotated with gene symbol. Accompanying tsv file of annotatated gene symbols in 
switch_vs_degs_RGC_vs
ROD209_vs_ROD45

`switch_RO_D209_vs_D45.pdf`: Switch plots of top genes for `RO_D209 vs D45` - genes that atleast one statistically significant switch. Consult top_genes.tsv for list of top 500 genes with switches.

`switch_RGC_vs_RO_D209.pdf`: Switch plots of top genes for `RGC vs RO_D209` - genes that atleast one statistically significant switch. Consult `top_genes.tsv` for list of top 500 genes with switches.


`switch_RGC_D209_vs_D45.pdf`: Switch plots of top genes for `RGC vs D45` - genes that atleast one statistically significant switch. Consult `top_genes.tsv` for list of top 500 genes with switches.


`isoform_volcano.pdf` : Pairwise Volcano plots of dIF, genes that have |dIF| > 0.5 and q value (FDR) < 0.05 are annotated with gene symbol. Accompanying tsvs of annotated genes -  `ROD209_vs_ROD45_switches.tsv`, `RGC_vs_ROD209_switches.tsv` , `RGC_vs_ROD45_switches.tsv`.  
