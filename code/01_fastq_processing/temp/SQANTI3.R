

output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/sqanti_output/EP1-BRN3B-RO_sqanti3_qc/"
sample <- "EP1-BRN3B-RO"
class.file <- paste0(output_dir, sample, "_", "classification.txt") 
report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
output_directory <- dirname(class.file)
output_name <- basename(report.prefix)
pdf.report.file <- paste0(report.prefix, "_SQANTI3_report.pdf");
class.file2 <- paste(report.prefix, "_classification_TPM.txt", sep='');



junc.file <- paste0(output_dir, sample, "_", "junctions.txt")
utilities.path <- "/users/sparthib/SQANTI3-5.2.1/utilities/"
saturation.curves <- FALSE
report.format <- "pdf"


#Packages 

library(ggplot2)
library(ggplotify)
library(scales)
library(reshape)
library(grid)
library(gridExtra)
library(NOISeq)
library(rmarkdown)
library(htmltools)
library(DT)
library(plyr)
library(plotly)
library(dplyr)

# PLOT length of isoforms
# p.length.all: length of all isoforms, regardless of category
# p.length.cat: length of isoforms, by category
# p.length.exon: length of isoforms, mono- vs mult-exon
# (optional) p.length.all.sample
# (optional) p.length.exon.sample

# p0: Splicing complexity (X) Isoforms per Gene (Y) Number of Genes
# p1: Distribution of Isoform Classification
# p2: Distribution of Ref Lengths (FSM ISM only)
# p3: Distribution of Ref Exons   (FSM ISM only)
# p4: Distribution of Isoform Lengths, by Classification
# p5: Distribution of Exon Counts, by Classification
# p6: Distribution of Mono- vs Multi-Exons, Novel vs Annotated Genes
# p7: Splicing complexity, Novel vs Annotated Genes
# p8: Transcript Expression (SR) by Structural Category { log2(TPM+1) }
# p9: Transcript Expression (FL) by Structural Category { log2(TPM+1) }
# p10: Gene Expression (SR), Annotated vs Novel genes { log2(Gene_TPM+1) }
# p11: Gene Expression (SR), Annotated vs Novel genes { log2(Gene_TPM+1) }
# p13: Gene Expression level in NNC/FSM containing genes
# p13.c: Gene Expression level in NNC/FSM containing genes

# p.classByLen.a: Structural categories with increasing transcript length, absolute
# p.classByLen.b: Structural categories with increasing transcript length, normalized
#
# p21.a: Distance to polyA site for FSM, absolute
# p21.b: Distance to polyA site for FSM, percentage
# p21.dist3.ISM.a: Distance to polyA site for ISM, absolute
# p21.dist3.ISM.b: Distance to polyA site for ISM, percentage
# p22.a: Distance to start site for FSM, absolute
# p22.b: Distance to start site for FSM, percentage

# p23.a: Splice Junctions by Classification (known canonical, known non-canonical, novel canonical, novel non-canonical)
# p23.b: Splice Junctions by Classification (canonical vs non-canonical)
#
# p28.a: Good Quality Control Attributes Across Structural Categories (annot, SJ, coverage)
# p28.aa: Good Quality Control Attributes Across Structural Categories (polyA, Cage)
# p28.a.SJ: Percentage of  All Canonical Junctions
# p28.a.Cov: Percentage of Splice Junctions With Short Reads Coverage
# p28.a.Cage: Percentage of Cage Support
# p28.a.polyA : Percentage of PolyA Support
# p28.a.annot : Percentage of Annotation Support
# p29.a: Splice Junction, % of RT switching, all junctions
# p29.b: Splice Junction, % of RT switching, unique junctions
#
# p30: intra-priming, by Classification
# p31: intra-priming, Mono- vs Multi-Exons
# p32: intra-priming, Coding vs Non-Coding


########## Classification information

data.class = read.table(class.file, header=T, as.is=T, sep="\t")
rownames(data.class) <- data.class$isoform

xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
subc.levels=c("alternative_3end",'alternative_3end5end', "alternative_5end","reference_match", "3prime_fragment","internal_fragment", "5prime_fragment","combination_of_known_junctions", "combination_of_known_splicesites", "intron_retention","no_combination_of_known_junctions", "mono-exon_by_intron_retention", "at_least_one_novel_splicesite", "mono-exon", "multi-exon")
subc.labels=c("Alternative 3'end", "Alternative 3'5'end", "Alternative 5'end", "Reference match", "3' fragment", "Internal fragment", "5' fragment", "Comb. of annot. junctions", "Comb. of annot. splice sites", "Intron retention", "Not comb. of annot. junctions", "Mono-exon by intron ret.", "At least 1 annot. don./accept.", "Mono-exon", "Multi-exon")
coding.levels=c("coding", "non_coding")
coding.labels=c("Coding", "Non coding")


data.class$structural_category = factor(data.class$structural_category,
                                        labels = xaxislabelsF1,
                                        levels = xaxislevelsF1,
                                        ordered=TRUE)
data.class$subcategory = factor(data.class$subcategory,
                                labels = subc.labels,
                                levels = subc.levels,
                                ordered=TRUE)
data.class$coding = factor(data.class$coding,
                           labels = coding.labels,
                           levels = coding.levels,
                           ordered=TRUE)
legendLabelF1 <- levels(as.factor(data.class$coding))


####### SRTM and SNTM functions

STM_function <- function(x){
  five=FALSE
  three=FALSE
  sj=FALSE
  ev_sj <- !is.na(x["min_cov"]) & as.numeric(x["exons"])>1
  ref_TSS <- FALSE
  ref_TTS <- FALSE
  
  if (!is.na(x["diff_to_gene_TSS"])){
    if (abs(as.numeric(x["diff_to_gene_TSS"]))<=50){
      ref_TSS <- TRUE
    }
  }
  
  if (!is.na(x["diff_to_gene_TTS"])){
    if (abs(as.numeric(x["diff_to_gene_TTS"]))<=50){
      ref_TTS <- TRUE
    }
  }
  
  w_cage <- !is.na(x["within_CAGE_peak"]) & x["within_CAGE_peak"]=="True"
  if ( ref_TSS | w_cage  ){
    five=TRUE
  }
  w_polya <- !is.na(x["within_polyA_site"]) & x["within_polyA_site"]=="True"
  if (ref_TTS | w_polya | !is.na(x["polyA_motif"])){
    three=TRUE
  }
  if (x["structural_category"]=="FSM" | x["structural_category"]=="ISM"){
    sj=TRUE
  }else if ( ev_sj ){
    if (as.numeric(x["min_cov"])>0){
      sj=TRUE
    }
  }

  if (five & three & sj){
    return("Fully supported")
  }else{
    return("Not fully supported")
  }
}

data.class$STM <- apply(data.class,1, STM_function)
################################

data.FSMISM <- subset(data.class, structural_category %in% c('FSM', 'ISM'))
data.NICNNC <- subset(data.class, structural_category %in% c("NIC", "NNC"))
data.other <- subset(data.class, structural_category %in% c("Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron"))
data.FSM <- subset(data.class, (structural_category=="FSM" & exons>1))
data.ISM <- subset(data.class, (structural_category=="ISM" & exons>1))
data.NNC <- subset(data.class, (structural_category=="NNC" & exons>1))
data.NIC <- subset(data.class, (structural_category=="NIC" & exons>1))
data.GenicGenomic <- subset(data.class, (structural_category=="Genic\nGenomic" & exons>1 ))
data.Antisense <- subset(data.class, (structural_category=="Antisense" & exons>1))
data.Fusion <- subset(data.class, (structural_category=="Fusion" & exons>1))
data.Intergenic <- subset(data.class, (structural_category=="Intergenic" & exons>1))
data.GenicIntron <- subset(data.class, (structural_category=="Genic\nIntron" & exons>1))

# subcategories data sets
#"FSM"
data.alt3end <- subset(data.FSM, (subcategory=="Alternative 3'end"))
data.alt35end <- subset(data.FSM, (subcategory=="Alternative 3'5'end"))
data.alt5end <- subset(data.FSM, (subcategory=="Alternative 5'end"))
data.refmatch <- subset(data.FSM, (subcategory=="Reference match"))
#"ISM"
data.3fragment <- subset(data.ISM, (subcategory=="3' fragment"))
data.int_fragment <- subset(data.ISM, (subcategory=="Internal fragment"))
data.5fragment <- subset(data.ISM, (subcategory=="5' fragment"))
data.intron_ret_ISM <- subset(data.ISM, (subcategory=="Intron retention"))
#"NIC"
data.comb_annot_js_NIC <- subset(data.NIC, (subcategory=="Comb. of annot. junctions"))
data.comb_annot_ss_NIC <- subset(data.NIC, (subcategory=="Comb. of annot. splice sites"))
data.intron_ret_NIC <- subset(data.NIC, (subcategory=="Intron retention"))
data.mono_ex_intron_ret_NIC <- subset(data.NIC, (subcategory=="Mono-exon by intron ret."))
#"NNC"
data.comb_annot_js_NNC <- subset(data.NNC, (subcategory=="Comb. of annot. junctions"))
data.comb_annot_ss_NNC <- subset(data.NNC, (subcategory=="Comb. of annot. splice sites"))
data.intron_ret_NNC <- subset(data.NNC, (subcategory=="Intron retention"))
data.mono_ex_intron_ret_NNC <- subset(data.NNC, (subcategory=="Mono-exon by intron ret."))
data.one_don_acc <- subset(data.NNC, (subcategory=="At least 1 annot. don./accept."))


subcategories.FSM <- list(data.alt3end, data.alt35end, data.alt5end, data.refmatch)
subcategories.ISM <- list(data.3fragment, data.int_fragment, data.5fragment, data.intron_ret_ISM)
subcategories.NIC <- list(data.comb_annot_js_NIC, data.comb_annot_ss_NIC, data.intron_ret_NIC, data.mono_ex_intron_ret_NIC)
subcategories.NNC <- list(data.comb_annot_js_NNC, data.comb_annot_ss_NNC, data.intron_ret_NNC, data.mono_ex_intron_ret_NNC, data.one_don_acc)

########### Junction information
data.junction <- read.table(junc.file, header=T, as.is=T, sep="\t")

# make a unique identifier using chrom_strand_start_end
data.junction$junctionLabel = with(data.junction, paste(chrom, strand, genomic_start_coord, genomic_end_coord, sep="_"))

data.junction$SJ_type <- with(data.junction, paste(junction_category,canonical,"SJ", sep="_"))
data.junction$SJ_type <- factor(data.junction$SJ_type, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T)

data.junction$structural_category = data.class[data.junction$isoform, "structural_category"]

uniqJunc <- unique(data.junction[,c("junctionLabel", "SJ_type", "total_coverage_unique")])
uniqJuncRTS <- unique(data.junction[,c("junctionLabel","SJ_type", "RTS_junction")])


########## Generating plots

#*** Global plot parameters

myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")
subcat.palette = c("Alternative 3'end"='#02314d',
                   "Alternative 3'5'end"='#0e5a87',
                   "Alternative 5'end"='#7ccdfc',
                   'Reference match'='#c4e1f2',
                   "3' fragment"='#c4531d',
                   "Internal fragment"='#e37744',  
                   "5' fragment"='#e0936e', 
                   "Comb. of annot. junctions"='#014d02',
                   "Comb. of annot. splice sites"='#379637',  
                   "Intron retention"='#81eb82', 
                   "Not comb. of annot. junctions"='#6ec091',
                   "Mono-exon by intron ret."='#4aaa72',
                   "At least 1 annot. don./accept."='#32734d',
                   "Mono-exon"='#cec2d2',
                   "Multi-exon"='#876a91')



cat.palette = c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")


mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=12) ) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size=12), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm"))


# Create a new attribute called "novelGene"

data.class$novelGene <- "Annotated Genes"
data.class[grep("novelGene", data.class$associated_gene), "novelGene"] <- "Novel Genes"
data.class$novelGene = factor(data.class$novelGene,
                              levels = c("Novel Genes","Annotated Genes"),
                              ordered=TRUE)

# Create a new attribute called "exonCat"

data.class[which(data.class$exons>1), "exonCat"] <- "Multi-Exon"
data.class[which(data.class$exons==1), "exonCat"] <- "Mono-Exon"
data.class$exonCat = factor(data.class$exonCat,
                            levels = c("Multi-Exon","Mono-Exon"),
                            ordered=TRUE)

canonical.labels=c("Canonical", "Non-canonical")
data.class$all_canonical = factor(data.class$all_canonical,
                                  labels=canonical.labels,
                                  levels = c("canonical","non_canonical"),
                                  ordered=TRUE)



if (!all(is.na(data.class$gene_exp))){
  isoPerGene = aggregate(data.class$isoform,
                         by = list("associatedGene" = data.class$associated_gene,
                                   "novelGene" = data.class$novelGene,
                                   "FSM_class" = data.class$FSM_class,
                                   "geneExp"=data.class$gene_exp),
                         length)
} else {
  isoPerGene = aggregate(data.class$isoform,
                         by = list("associatedGene" = data.class$associated_gene,
                                   "novelGene" = data.class$novelGene,
                                   "FSM_class" = data.class$FSM_class),
                         length)
}
# assign the last column with the colname "nIso" (number of isoforms)
colnames(isoPerGene)[ncol(isoPerGene)] <- "nIso"


isoPerGene$FSM_class2 = factor(isoPerGene$FSM_class,
                               levels = c("A", "B", "C"),
                               labels = c("MonoIsoform Gene", "MultiIsoform Genes\nwithout expression\nof a FSM", "MultiIsoform Genes\nexpressing at least\none FSM"),
                               ordered=TRUE)

isoPerGene$novelGene = factor(isoPerGene$novelGene,
                              levels = c("Annotated Genes", "Novel Genes"),
                              ordered=TRUE)

max_iso_per_gene <- max(isoPerGene$nIso)
if (max_iso_per_gene >= 6) {
  isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3,5,max_iso_per_gene+1), labels = c("1", "2-3", "4-5", ">=6"));
} else if (max_iso_per_gene >= 5) {
  isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3,5), labels = c("1", "2-3", "4-5"));
} else if (max_iso_per_gene >= 3) {
  isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3), labels = c("1", "2-3"));
} else {
  isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1), labels = c("1"));
}

########################################
######### LENGTH PLOTS  ################
########################################

p.length.all <- ggplot(data.class, aes(x=length)) +
  geom_histogram(binwidth=100) +
  labs(x="Transcript length", y="Count", title="All Transcript Lengths Distribution") +
  theme(legend.position="top") +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  mytheme

p.length.cat <- ggplot(data.class, aes(x=length, color=structural_category)) +
  geom_freqpoly(binwidth=100, size=1) +
  scale_color_manual(values = cat.palette, name="Structural Category") +
  labs(x="Transcript length", y="Count", title="Transcript Lengths Distribution by Structural Category") +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  mytheme+
  theme(legend.position="bottom", legend.title=element_blank())

p.length.exon <- ggplot(data.class, aes(x=length, color=exonCat)) +
  geom_freqpoly(binwidth=100, size=1) +
  labs(x="Transcript length", y="Count", title="Mono- vs Multi- Exon Transcript Lengths Distribution") +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  mytheme +
  theme(legend.position="bottom", legend.title=element_blank()) 


#**** PLOT 0: Distribution of Number of Isoforms per Gene

p0 <- ggplot(isoPerGene, aes(x=nIsoCat, fill=nIsoCat)) +
  geom_bar(stat="count", aes(y= (..count..)/sum(..count..)*100), color="black", size=0.3, width=0.7) +
  guides(fill="none") +
  scale_y_continuous(labels = function(x) format(x), expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(2:5)]) +
  labs(x ="Isoforms per gene", title="Number of Isoforms per Gene\n\n\n", y = "Genes, %") +
  mytheme

#**** PLOT 1: Structural Classification

p1 <- ggplot(data=data.class, aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100, alpha=coding, fill=structural_category), color="black", size=0.3, width=0.7) +
  #geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25)  +
  scale_x_discrete(drop=FALSE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Coding prediction",
                     labels = legendLabelF1)+
  xlab("") +
  ylab("Transcripts, %") +
  mytheme +
  geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))

p1.s.titles = list("Isoform Distribution Across FSM\n\n",
                   "Isoform Distribution Across ISM\n\n",
                   "Isoform Distribution Across NNC\n\n",
                   "Isoform Distribution Across NIC\n\n",
                   "Isoform Distribution Across Genic Genomic\n\n",
                   "Isoform Distribution Across Antisense\n\n",
                   "Isoform Distribution Across Fusion\n\n",
                   "Isoform Distribution Across Intergenic\n\n",
                   "Isoform Distribution Across Genic Intron\n\n")

categories.list=list(data.FSM, data.ISM, data.NNC, data.NIC, data.GenicGenomic, data.Antisense, 
                     data.Fusion, data.Intergenic, data.GenicIntron)

p1.s.list = list()
for(i in 1:length(categories.list)){
  c<-data.frame(categories.list[i])
  if (!(dim(c))[1]==0){
    p1.s <- ggplot(data=c, aes(x=subcategory)) +
      geom_bar(aes(y = (..count..)/sum(..count..)*100, alpha=coding, fill=subcategory), color="black", size=0.3, width=0.7) +
      scale_x_discrete(drop=TRUE) +
      scale_alpha_manual(values=c(1,0.3), name = "Coding prediction", labels = legendLabelF1)+
      ylab("Transcripts, %") +
      mytheme +
      geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values = subcat.palette, guide='none') +
      ggtitle(p1.s.titles[i])+
      scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
      theme(axis.title.x=element_blank())  
    p1.s.list[[i]] = p1.s
  }
}


#**** PLOTS 2-3: refLength and refExons for ISM and FSM transcripts. Plot if any ISM or FSM transcript

if (nrow(data.FSMISM) > 0) {
  
  p2 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=ref_length/1000, fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
    scale_fill_manual(values = myPalette) +
    scale_x_discrete(drop = TRUE) +
    guides(fill="none") +
    xlab("") +
    ylab("Matched reference length, kb") +
    labs(title="Length Distribution of Matched Reference Transcripts\n\n\n",
         subtitle="Applicable Only to FSM and ISM Categories\n\n")
  
  p3 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=ref_exons, fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    scale_x_discrete(drop = TRUE) +
    xlab("") +
    ylab("Matched reference exon count") +
    scale_fill_manual(values = myPalette) +
    guides(fill="none") +
    mytheme +
    labs(title="Exon Count Distribution of Matched Reference Transcripts\n\n\n",
         subtitle="Applicable Only to FSM and ISM Categories\n\n")
  
}


#****  PLOT 4: Transcript lengths by category

p4 <- ggplot(data=data.class, aes(x=structural_category, y=length, fill=structural_category)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=FALSE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = cat.palette) +
  guides(fill="none") +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Transcript Lengths by Structural Classification\n\n" ) +
  theme(axis.title.x=element_blank())

p4.s1 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=length, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=TRUE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = subcat.palette, drop=TRUE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position="right", legend.title=element_blank()) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Transcript Lengths by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())

p4.s2 <- ggplot(data=data.NICNNC, aes(x=structural_category, y=length, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=TRUE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = subcat.palette, drop=TRUE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position="right", legend.title=element_blank())+
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Transcript Lengths by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())

p4.s3 <- ggplot(data=data.other, aes(x=structural_category, y=length, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=TRUE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = subcat.palette, drop=TRUE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position="right", legend.title=element_blank())+
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Transcript Lengths by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())


##**** PLOT 5: Exon counts by category

p5 <- ggplot(data=data.class, aes(x=structural_category, y=exons, fill=structural_category)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Number of exons") +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(values = cat.palette) +
  guides(fill="none") +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Exon Counts by Structural Classification\n\n" ) +
  theme(axis.title.x=element_blank())

###Exon counts by subcategory
p5.s1 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=exons, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  theme(legend.position="right", legend.title=element_blank()) +
  theme(axis.title.x=element_blank())+
  ylab("Number of exons") +
  scale_x_discrete(drop=TRUE) +
  scale_fill_manual(values = subcat.palette) +
  ggtitle("Exon Counts by Subcategory\n\n" )


p5.s2 <- ggplot(data=data.NICNNC, aes(x=structural_category, y=exons, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Number of exons") +
  scale_x_discrete(drop=TRUE) +
  scale_fill_manual(values = subcat.palette) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  theme(legend.position="right", legend.title=element_blank())+
  ggtitle("Exon Counts by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())

p5.s3 <- ggplot(data=data.other, aes(x=structural_category, y=exons, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Number of exons") +
  scale_x_discrete(drop=TRUE) +
  scale_fill_manual(values = subcat.palette) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  theme(legend.position="right", legend.title=element_blank())+
  ggtitle("Exon Counts by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())


##### STM plots

data.FSMISMNICNNC <- rbind(data.FSMISM, data.NICNNC )
pSTM <- ggplot(data=data.FSMISMNICNNC, aes(x=structural_category)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=structural_category), color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, count") +
  mytheme +
  geom_blank(aes(y=(..count..)), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  theme(legend.position = "right")


pSTM_perc <- ggplot(data=data.FSMISMNICNNC, aes(x=structural_category)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=structural_category), position="fill", color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, %") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1)), labels = scales::percent) +
  theme(legend.position = "right")


pSTM.s1 <- ggplot(data=data.FSMISM, aes(x=subcategory)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, count") +
  mytheme +
  facet_grid(.~ structural_category, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = subcat.palette, guide = "none") +
  ggtitle("Isoform Distribution Across Structural Subcategories\n\n",
          subtitle = "FSM and ISM" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x=element_text(size=10)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  theme(legend.position = "right")

pSTM_perc.s1 <- ggplot(data=data.FSMISM, aes(x=subcategory)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), position="fill", color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, %") +
  mytheme +
  facet_grid(.~ structural_category, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = subcat.palette, guide = "none") +
  ggtitle("Isoform Distribution Across Structural Subcategories\n\n",
          subtitle = "FSM and ISM" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(size=10)) +
  scale_y_continuous(expand=expansion(mult = c(0,0)), labels = scales::percent) +
  theme(legend.position = "right")


pSTM.s2 <- ggplot(data=data.NICNNC, aes(x=subcategory)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, count") +
  mytheme +
  facet_grid(.~ structural_category, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = subcat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Subcategories\n\n",
          subtitle = "NIC and NNC") +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(size=10)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  theme(legend.position = "right")

pSTM_perc.s2 <- ggplot(data=data.NICNNC, aes(x=subcategory)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), position="fill", color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, %") +
  mytheme +
  facet_grid(.~ structural_category, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = subcat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Subcategories\n\n",
          subtitle = "NIC and NNC" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(size=10)) +
  scale_y_continuous(expand=expansion(mult = c(0,0)), labels = scales::percent) +
  theme(legend.position = "right")




##**** PLOT 6: Mono vs Multi-exon distribution for Known vs Novel Genes

p6 <- ggplot(data=data.class, aes(x=novelGene)) +
  geom_bar(position="fill",aes(y = (..count..)/sum(..count..), fill=exonCat), color="black", size=0.3, width=0.5) +
  scale_x_discrete(drop=FALSE) +
  scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                     labels=c("0","25","50","75","100"), expand = c(0,0)) +
  scale_fill_manual(name = "Transcript type",
                    values = myPalette[c(2:5)]) +
  ylab("Transcripts, %") +
  mytheme +
  theme(axis.title.x=element_blank()) +
  theme(legend.position="bottom") +
  ggtitle("Distribution of Mono- vs Multi-Exon Transcripts\n\n" )

# p7: Distribution of Number of Isoforms, separated by Novel vs Annotated Genes

p7 <- ggplot(data=isoPerGene, aes(x=novelGene)) +
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=nIsoCat), color="black", size=0.3, width=0.5) +
  scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                     labels=c("0","25","50","75","100"), expand = c(0,0)) +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(name = "Isoforms Per Gene",
                    values = myPalette[c(2:5)]) +
  ylab("Genes, %") +
  xlab("Gene Type") +
  mytheme +
  theme(axis.title.x=element_blank()) +
  theme(legend.position="bottom") +
  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
  labs(title="Number of Isoforms per Gene\n\n\n",
       subtitle="Known vs Novel Genes\n\n")


##**** PLOT  absolute and normalized % of different categories with increasing transcript length
# requires use of dplyr package
data.class$lenCat <- as.factor(as.integer(data.class$length %/% 1000))
data.class.byLen <- data.class %>% dplyr::group_by(lenCat, structural_category) %>% dplyr::summarise(count=dplyr::n() ) %>% mutate(perc=count/sum(count))
data.class.byLen$structural_category <- factor(data.class.byLen$structural_category, levels=(xaxislabelsF1), order=TRUE)

p.classByLen.a <- ggplot(data.class.byLen, aes(x=lenCat, y=count, fill=factor(structural_category))) +
  geom_bar(stat='identity', color="black", size=0.3, width=0.85) +
  scale_fill_manual(values = cat.palette, guide='none', name="Structural Category") +
  mytheme+
  theme(legend.position="right") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
  labs(x="Transcript length, kb", y="Counts", title="Structural Categories by Transcript Length")


p.classByLen.b <- ggplot(data.class.byLen, aes(x=lenCat, y=perc*100, fill=factor(structural_category))) +
  geom_bar(stat='identity', color ="black", size=0.3, width=0.85) +
  scale_fill_manual(values = cat.palette, guide='none', name="Structural Category") +
  mytheme+
  theme(legend.position="right")  +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
  labs(x="Transcript length, kb", y="%", title="Structural Categories by Transcript Length\n\n\n")

##**** PLOT 8: Expression, if isoform expression provided (iso_exp is in TPM)
if (!all(is.na(data.class$iso_exp))){
  p8 <- ggplot(data=data.class, aes(x=structural_category, y=log2(iso_exp+1), fill=structural_category)) +
    geom_boxplot(color="black", size=0.3,  outlier.size = 0.2) +
    scale_x_discrete(drop=FALSE) +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = cat.palette) +
    guides(fill="none") +
    mytheme  + theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Transcript Expression by Structural Category\n\n" )
}
###Expression, if isoform expression provided (iso_exp is in TPM) by subcategory
if (!all(is.na(data.FSMISM$iso_exp))){
  p8.s1 <- ggplot(data=data.FSMISM, aes(x=subcategory, y=log2(iso_exp+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3,  outlier.size = 0.2, position="dodge") +
    scale_x_discrete(drop=TRUE) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme  + theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Transcript Expression by Subcategory\n\n" )
}

if (!all(is.na(data.NICNNC$iso_exp))){
  p8.s2 <- ggplot(data=data.NICNNC, aes(x=subcategory, y=log2(iso_exp+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3,  outlier.size = 0.2, position="dodge") +
    scale_x_discrete(drop=TRUE) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme  + theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Transcript Expression by Subcategory\n\n" )
}

if (!all(is.na(data.other$iso_exp))){
  p8.s3 <- ggplot(data=data.other, aes(x=subcategory, y=log2(iso_exp+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3,  outlier.size = 0.2, position="dodge") +
    scale_x_discrete(drop=TRUE) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme  + theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Transcript Expression by Subcategory\n\n" )
}


# PLOT 9: FL number, if FL count provided
# convert FL count to TPM
if (!all(is.na(data.class$FL))){
  p9 <- ggplot(data=data.class, aes(x=structural_category, y=log2(FL_TPM+1), fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size=0.2) +
    ylab("log2(FL_TPM+1)") +
    scale_x_discrete(drop=FALSE) +
    scale_fill_manual(values = cat.palette) +
    guides(fill="none") +
    mytheme +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Long Reads Count by Structural Category\n\n" )
}

if (!all(is.na(data.FSMISM$FL))){
  p9.s1 <- ggplot(data=data.FSMISM, aes(x=subcategory, y=log2(FL_TPM+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3, outlier.size=0.1) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(FL_TPM+1)") +
    scale_x_discrete(drop=TRUE) +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Long Reads Count by Subcategory\n\n" )
}

if (!all(is.na(data.NICNNC$FL))){
  p9.s2 <- ggplot(data=data.NICNNC, aes(x=subcategory, y=log2(FL_TPM+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3, outlier.size=0.1) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(FL_TPM+1)") +
    scale_x_discrete(drop=TRUE) +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Long Reads Count by Subcategory\n\n" )
}

if (!all(is.na(data.other$FL))){
  p9.s3 <- ggplot(data=data.other, aes(x=subcategory, y=log2(FL_TPM+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3, outlier.size=0.1) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(FL_TPM+1)") +
    scale_x_discrete(drop=TRUE) +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Long Reads Count by Subcategory\n\n" )
}



# PLOT 10: Gene Expression, if expresion provided
if (!all(is.na(data.class$iso_exp))){
  p10 <- ggplot(data=isoPerGene, aes(x=novelGene, y=log2(geneExp+1), fill=novelGene)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    scale_x_discrete(drop=FALSE) +
    xlab("Structural classification") +
    ylab("log2(Gene_TPM+1)") +
    scale_fill_manual(values = myPalette[c(3:4)]) +
    guides(fill="none") +
    mytheme +
    theme(axis.title.x=element_blank()) +
    ggtitle("Annotated vs Novel Gene Expression\n\n" )
}


# PLOT 11: Gene FL number, if FL count provided

if (!all(is.na(data.class$FL))){
  FL_gene <- aggregate(as.integer(data.class$FL), by = list("associatedGene" = data.class$associated_gene), sum)
  colnames(FL_gene)[ncol(FL_gene)] <- "FL_gene"
  isoPerGene <- merge(isoPerGene, FL_gene, by="associatedGene")
  total_fl <- sum(data.class$FL, na.rm=T)
  isoPerGene$FL_gene_TPM <- isoPerGene$FL_gene*(10**6)/total_fl
  
  p11 <- ggplot(data=isoPerGene, aes(x=novelGene, y=log2(FL_gene_TPM+1), fill=novelGene)) +
    geom_boxplot(color="black", size=0.3,outlier.size = 0.2) +
    scale_x_discrete(drop=FALSE) +
    ylab("log2(FL_TPM+1)") +
    scale_fill_manual(values = myPalette[c(3:4)]) +
    guides(fill="none") +
    mytheme +
    theme(axis.title.x=element_blank()) +
    ggtitle("Number of FL reads per Gene by Type of Gene Annotation\n\n" )
  
}




# PLOT 12: NNC expression genes vs not NNC expression genes
# NNC expression genes vs not NNC expression genes

if (!all(is.na(data.class$gene_exp))){
  if (nrow(data.class[data.class$structural_category=="NNC",])!=0){
    
    NNC_genes <- unique(data.class[data.class$structural_category=="NNC","associated_gene"])
    notNNC_genes <- unique(data.class[!data.class$associated_gene%in%NNC_genes,"associated_gene"])
    isoPerGene[isoPerGene$associatedGene %in% notNNC_genes, "NNC_class"] <- "Genes without\n NNC isoforms"
    isoPerGene[isoPerGene$associatedGene %in% NNC_genes, "NNC_class"] <- "Genes with\n NNC isoforms"
    
    isoPerGene$NNC_class <- factor(isoPerGene$NNC_class, levels=c("Genes with\n NNC isoforms","Genes without\n NNC isoforms"),
                                   labels=c("Genes with\n NNC isoforms","Genes without\n NNC isoforms"), order=T)
    
    p12 <- ggplot(data=isoPerGene[!is.na(isoPerGene$NNC_class),], aes(x=NNC_class, y=log2(geneExp+1), fill=NNC_class)) +
      geom_boxplot(color="black", size=0.3, outlier.size=0.2) +
      xlab("") +  
      ylab("log2(Gene_TPM+1)") +
      scale_x_discrete(drop=FALSE) +
      scale_fill_manual(values = c(myPalette[4],"grey38")) +
      guides(fill="none") +
      mytheme +
      theme(axis.title.x=element_blank()) + 
      ggtitle("Gene Expression of NNC And Not NNC Containing Genes\n\n" )
  }
}


if (!all(is.na(data.class$gene_exp))){
  if (nrow(data.class[data.class$structural_category=="NNC",])!=0 & nrow(data.class[data.class$structural_category=="FSM",])!=0 ){
    
    FSM_just_genes = unique(data.class[data.class$FSM_class=="A" & data.class$structural_category=="FSM","associated_gene"])
    NNC_just_genes = unique(data.class[data.class$FSM_class=="A" & data.class$structural_category=="NNC","associated_gene"])
    FSMandNNCgenes = unique(data.class[data.class$FSM_class=="C" & data.class$structural_category=="NNC","associated_gene"])
    isoPerGene[isoPerGene$associatedGene %in% FSMandNNCgenes, "FSM_NNC_class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
    isoPerGene[isoPerGene$associatedGene %in% NNC_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only NNC isoforms"
    isoPerGene[isoPerGene$associatedGene %in% FSM_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only FSM isoforms"
    data.class[data.class$associated_gene %in% FSMandNNCgenes, "class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
    data.class[data.class$associated_gene %in% NNC_just_genes, "class"] <- "Genes expressing\n only NNC isoforms"
    data.class[data.class$associated_gene %in% FSM_just_genes, "class"] <- "Genes expressing\n only FSM isoforms"
    
    isoPerGene$FSM_NNC_class = factor(isoPerGene$FSM_NNC_class, levels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"),
                                      labels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"), order=T)
    
    p13 <- ggplot(data=isoPerGene[!is.na(isoPerGene$FSM_NNC_class),], aes(x=FSM_NNC_class, y=log2(geneExp+1), fill=FSM_NNC_class)) +
      geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
      ylab("log2( # Short reads per gene + 1)") +
      theme(axis.title.x=element_blank()) +
      #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
      scale_fill_manual(values = c("grey38",myPalette[[4]],myPalette[[1]])) +
      guides(fill="none") +
      mytheme +
      theme(axis.title.x=element_blank()) +
      ggtitle("Gene Expression Level in NNC/FSM Containing Genes\n\n" ) +
      scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                                "Genes expressing\n only FSM isoforms",
                                "Genes expressing\n only NNC isoforms"),
                       labels=c("NNC/FSM genes",
                                "FSM genes",
                                "NNC genes"), drop=FALSE) 
    p13.c <- ggplot(data=data.class[!is.na(data.class$class),], aes(x=class, y=log2(iso_exp+1), fill=structural_category)) +
      geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
      ylab("log2( # Short reads per transcript + 1)") +
      theme(axis.title.x=element_blank()) +
      scale_fill_manual(values = myPalette) +
      guides(fill="none") +
      mytheme +
      theme(axis.title.x=element_blank()) +
      ggtitle("Transcript Expression Level in NNC/FSM Containing Genes\n\n" ) +
      scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                                "Genes expressing\n only FSM isoforms",
                                "Genes expressing\n only NNC isoforms"),
                       labels=c("NNC/FSM genes",
                                "FSM genes",
                                "NNC genes"), drop=F) 
    
  }
}
    

# PLOT 23: Junction categories


if (nrow(data.junction) > 0){
  data.junction$junctionLabel = with(data.junction, paste(chrom, strand,genomic_start_coord, genomic_end_coord, sep="_"))
  
  data.junction$canonical_known = with(data.junction, paste(junction_category,canonical,"SJ", sep="_"))
  data.junction$canonical_known=as.factor(data.junction$canonical_known)
  data.junction$canonical_known = factor(data.junction$canonical_known, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                         labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T) 
  data.junction$structural_category = data.class[data.junction$isoform,"structural_category"]
  ##    data.junction$TSSrange =cut(data.junction$transcript_coord, breaks = c(0, 40, 80, 120, 160, 200, 10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))
  
  p23.a <- ggplot(data.junction, aes(x=structural_category)) +
    geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=SJ_type), color="black",  size=0.3, width = 0.7) +
    scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                       labels=c("0","25","50","75","100"), expand = c(0,0)) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
    ylab("Splice junctions, %") +
    mytheme +
    guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
    theme(legend.position="bottom", legend.title=element_blank())  +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Distribution of Splice Junctions by Structural Classification\n\n\n")
  
  t <- subset(data.class, exons > 1)  # select only multi-exon isoforms
  
  p23.b <- ggplot(data=t, aes(x=structural_category)) +
    geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=all_canonical), color="black", size=0.3, width = 0.7) +
    scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                       labels=c("0","25","50","75","100"), expand = c(0,0)) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
    xlab("") +
    ylab("Transcripts, %") +
    mytheme +
    guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
    theme(legend.position="bottom", legend.title=element_blank())  +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Distribution of Transcripts by Splice Junctions\n\n\n")
}
  
  
### Bad quality control attributes
if (nrow(data.junction) > 0){
  t3.data.sets <- list()
  t3.list <- list()
  # (Fran) ToDo: USE COVERAGE DATA LATER
  # for FSM, ISM, NIC, and NNC, plot the percentage of RTS and non-canonical junction
  x <- filter(data.class, structural_category %in% c("FSM", "ISM", "NIC", "NNC" ) & exons > 1)
  
  t1.RTS <- group_by(x, structural_category, RTS_stage) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  t2.RTS <- group_by(x, structural_category) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  t3.RTS <- merge(t1.RTS, t2.RTS, by="structural_category")
  t3.RTS <- t3.RTS[-which(t3.RTS$structural_category=="ISM"),]
  t3.RTS$perc <- t3.RTS$count.x / t3.RTS$count.y * 100
  t3.RTS <- subset(t3.RTS, RTS_stage=='TRUE');
  n_t3.RTS <- dim(t3.RTS)[1];
  if (n_t3.RTS > 0) {
    t3.RTS$Var <- "RT switching"
  }
  
  # Liz: this is a placeholder for dealing with all_canonical being NA instead of "Non-canonical"
  x[is.na(x$all_canonical), "all_canonical"] <- "Non-canonical"
  t1.SJ <- group_by(x, structural_category, all_canonical) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  t3.SJ <- merge(t1.SJ, t2.RTS, by="structural_category")
  t3.SJ$perc <- t3.SJ$count.x / t3.SJ$count.y * 100
  t3.a.SJ <- subset(t3.SJ, all_canonical=='Canonical');
  t3.SJ <- subset(t3.SJ, all_canonical=='Non-canonical');
  n_t3.SJ <- dim(t3.SJ)[1];
  if (n_t3.SJ > 0) {
    t3.SJ$Var <- "Non-canonical"
    t3.a.SJ$Var <- 'Canonical'
  }
  
  if (!all(is.na(x$predicted_NMD))){
    x[which(x$predicted_NMD=="TRUE"),"predicted_NMD"]="Predicted NMD"
    x[which(x$predicted_NMD=="FALSE"),"predicted_NMD"]="Not NMD predicted"
    t1.NMD <- group_by(x, structural_category, predicted_NMD) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
    t3.NMD <- merge(t1.NMD, t2.RTS, by="structural_category")
    t3.NMD$perc <- t3.NMD$count.x / t3.NMD$count.y * 100
    t3.NMD <- subset(t3.NMD, predicted_NMD=='Predicted NMD');
    t3.NMD$Var=t3.NMD$predicted_NMD
  }
  if (!all(is.na(x$min_cov))){
    x[which(x$min_cov==0),"Coverage_SJ"]="Not Coverage SJ"
    x[which(x$min_cov>0),"Coverage_SJ"]="Coverage SJ"
    t1.Cov <- group_by(x, structural_category, Coverage_SJ) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
    t3.Cov <- merge(t1.Cov, t2.RTS, by="structural_category")
    t3.Cov$perc <- t3.Cov$count.x / t3.Cov$count.y * 100
    t3.a.Cov <- subset(t3.Cov, Coverage_SJ=='Coverage SJ');
    t3.Cov <- subset(t3.Cov, Coverage_SJ=='Not Coverage SJ');
    t3.Cov$Var=t3.Cov$Coverage_SJ
    t3.a.Cov$Var=t3.a.Cov$Coverage_SJ
    t3.data.sets[[length(t3.data.sets) + 1]]=x$min_cov
    t3.list[[length(t3.list) + 1]]=t3.a.Cov
  }
  
  # Check for Non-empty Data
  if (nrow(x) > 0) {
    x[which(x$diff_to_gene_TSS<=50),"Annotation"] <- "Annotated"
    x[which(x$diff_to_gene_TSS>50),"Annotation"] <- "Not annotated"
    t1.annot <- group_by(x, structural_category, Annotation) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
    t3.annot <- merge(t1.annot, t2.RTS, by="structural_category")
    t3.annot$perc <- t3.annot$count.x / t3.annot$count.y * 100
    t3.annot <- subset(t3.annot, Annotation=='Annotated');
    t3.annot$Var=t3.annot$Annotation
    p28.a.annot <- ggplot(t3.annot, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[6] ,color="black") +
      geom_text(label=paste(round(t3.annot$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Annotation Support\n\n") 
  }
  
  p28.RTS <- ggplot(t3.RTS, aes(x=structural_category, y=perc)) +
    geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[11], color="black") +
    geom_text(label=paste(round(t3.RTS$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
    scale_fill_manual(values = myPalette[9:11]) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    ylab("Isoforms, %") +
    xlab("") +
    mytheme +
    theme(legend.position="bottom", axis.title.x = element_blank()) +
    ggtitle("RT-switching\n\n")
  
  p28.SJ <- ggplot(t3.SJ, aes(x=structural_category, y=perc)) +
    geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[9] ,color="black") +
    geom_text(label=paste(round(t3.SJ$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
    scale_fill_manual(values = myPalette[9:11]) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    ylab("Isoforms, %") +
    xlab("") +
    mytheme +
    theme(legend.position="bottom", axis.title.x = element_blank()) +
    ggtitle("Non-Canonical Junctions\n\n")
  p28.a.SJ <- ggplot(t3.a.SJ, aes(x=structural_category, y=perc)) +
    geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[7] ,color="black") +
    geom_text(label=paste(round(t3.a.SJ$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    ylab("Isoforms, %") +
    xlab("") +
    mytheme +
    theme(legend.position="bottom", axis.title.x = element_blank()) +
    ggtitle("All Canonical Junctions\n\n")
  
  if (n_t3.SJ>0 & n_t3.RTS>0 & !all(is.na(x$min_cov)) & all(is.na(x$predicted_NMD))){
    p28.Cov <- ggplot(t3.Cov, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
      geom_text(label=paste(round(t3.Cov$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Splice Junctions Without Short Reads Coverage\n\n")
    p28.a.Cov <- ggplot(t3.a.Cov, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
      geom_text(label=paste(round(t3.a.Cov$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Splice Junctions With Short Reads Coverage\n\n") 
    
    
    t3 <- rbind(t3.RTS[,c(1,5,6)],t3.SJ[,c(1,5,6)], t3.Cov[,c(1,5,6)])
    
    p28 <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Transcripts, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Summary Features of Bad Quality\n\n" ) +
      theme(legend.title = element_blank())
    
    #good quality control
    t3.a <- rbind(t3.annot[,c(1,5,6)], t3.a.SJ[,c(1,5,6)], t3.a.Cov[,c(1,5,6)])
    
  }else if (n_t3.SJ>0 & n_t3.RTS>0 & all(is.na(x$min_cov)) & all(is.na(x$predicted_NMD))) {
    t3=rbind(t3.RTS[,c(1,5,6)],t3.SJ[,c(1,5,6)])
    p28 <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[c(9,11)]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Transcripts, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Quality Control Attributes Across Structural Categories\n\n" ) +
      theme(legend.title = element_blank())
    #good quality control
    t3.a=rbind(t3.annot[,c(1,5,6)], t3.a.SJ[,c(1,5,6)])
    
  }else if (n_t3.SJ>0 & n_t3.RTS>0 & all(is.na(x$min_cov)) & !all(is.na(x$predicted_NMD))){
    p28.NMD <- ggplot(t3.NMD, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[5], color="black") +
      geom_text(label=paste(round(t3.NMD$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Nonsense-Mediated Decay by Structural Category\n\n")
    t3=rbind(t3.RTS[,c(1,5,6)],t3.SJ[,c(1,5,6)], t3.NMD[,c(1,5,6)])
    p28 <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[c(9,5,11)]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Transcripts, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Quality Control Attributes Across Structural Categories\n\n" ) +
      theme(legend.title = element_blank())
    #good quality control
    t3.a=rbind(t3.annot[,c(1,5,6)], t3.a.SJ[,c(1,5,6)])
    
  }else if (n_t3.SJ>0 & n_t3.RTS>0) {
    p28.NMD <- ggplot(t3.NMD, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[5], color="black") +
      geom_text(label=paste(round(t3.NMD$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Nonsense-Mediated Decay by Structural Category\n\n")
    p28.Cov <- ggplot(t3.Cov, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
      geom_text(label=paste(round(t3.Cov$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Splice Junctions Without Short Read Coverage\n\n")
    p28.a.Cov <- ggplot(t3.a.Cov, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
      geom_text(label=paste(round(t3.a.Cov$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Splice Junctions With Short Read Coverage\n\n")
    t3=rbind(t3.RTS[,c(1,5,6)],t3.SJ[,c(1,5,6)],t3.Cov[,c(1,5,6)], t3.NMD[,c(1,5,6)])
    
    p28 <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[c(9,10,5,11)]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Transcripts, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Quality Control Attributes Across Structural Categories\n\n" ) +
      theme(legend.title = element_blank())
    #good quality control
    t3.a <- rbind(t3.annot[,c(1,5,6)], t3.a.SJ[,c(1,5,6)], t3.a.Cov[,c(1,5,6)])
    
  }
  
}

# Check if t3.SJ is not empty
if (nrow(t3.SJ) > 0) {
  t3.aa <-  rbind(t3.annot[,c("structural_category", "perc", "Var")], t3.a.SJ[,c(1,5,6)])
  
  for(i in 1:length(t3.list)){
    set=data.frame(t3.data.sets[i])
    c=data.frame(t3.list[i])
    if (!all(is.na(set))){
      t.temp=t3.aa
      t3.aa = rbind(t.temp, c[,c(1,5,6)])
    }
  }
  
  p28.a <- ggplot(data=t3.aa, aes(x=structural_category, y=perc, fill= Var)) +
    geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    scale_fill_manual(values = c(myPalette)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    ylab("Transcripts, %") +
    xlab("") +
    mytheme +
    theme(legend.position="bottom", axis.title.x = element_blank()) +
    ggtitle( "Good Quality Control Attributes Across Structural Categories\n\n" ) +
    theme(axis.text.y = element_text(size=10),
          axis.text.x  = element_text(size=10))+
    theme(legend.title = element_blank())
}
    


### GENERATE PDF REPORT ###

generatePDFreport = function() 
{
  pdf(file=pdf.report.file, width = 7.5, height = 6.5)
  
  #cover
  grid.newpage()
  cover <- textGrob("SQANTI3 report",
                    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
  grid.draw(cover)
  
  # TABLE 1: Number of isoforms in each structural category
  
  freqCat <- as.data.frame(table(data.class$structural_category))
  #freqCat$ranking = order(freqCat$Freq,decreasing = T)
  table1 <- tableGrob(freqCat, rows = NULL, cols = c("Category","Isoforms, count"))
  title1 <- textGrob("Transcript Classification\n", gp=gpar(fontface="italic", fontsize=17), vjust = -
                       3.2)
  gt1 <- gTree(children=gList(table1, title1))
  
  
  # TABLE 2: Number of Novel vs Known Genes
  freqCat = as.data.frame(table(isoPerGene$novelGene))
  table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","Genes, count"))
  title2 <- textGrob("Gene Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
  gt2 <- gTree(children=gList(table2, title2))
  
  # TABLE 3: Junction Classification
  
  uniq_sj_count <- nrow(uniqJunc)
  
  freqCat <- as.data.frame(table(uniqJunc$SJ_type))
  freqCat$Var1 <- gsub(" ", "", freqCat$Var1)
  freqCat$Var1 <- gsub("\n", " ", freqCat$Var1)
  freqCat$Frac <- round(freqCat$Freq*100 / uniq_sj_count, 2)
  table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","SJs, count","Percent"))
  title2 <- textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
  gt3 <- gTree(children=gList(table2, title2))
  
  
  # TABLE 4: Summary number of Unique Isoforms and Unique Genes
  nGenes = nrow(isoPerGene)
  nIso = nrow(data.class)
  sn = paste("Unique Genes: ", nGenes, "\n", "Unique Isoforms: ", nIso)
  gt4 <- textGrob(sn, gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  
  
  # Plot Table 1 and Table 2
  grid.arrange(gt4,gt2,gt1, layout_matrix = cbind(c(1,2),c(1,4)))
  grid.arrange(gt3)
  
  s <- textGrob("Gene Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p0)
  print(p7)
  print(p6)
  print(p.classByLen.a)
  print(p.classByLen.b)
  
  if (!all(is.na(data.class$iso_exp))){
    print(p10)
  }
  if (!all(is.na(data.class$FL))){
    print(p11)
  }
  
  # PLOT length of isoforms
  # p.length.all: length of all isoforms, regardless of category
  # p.length.cat: length of isoforms, by category
  # p.length.exon: length of isoforms, mono- vs mult-exon/ufrc/conesa/fpardopalacios/SQANTI_QDE/SQANTI3/melanoma_example/melanoma_chr13_tappAS_annot_from_SQANTI3.gff3
  # (optional) p.length.all.sample: length of all isoforms by sample
  print(p.length.all)
  print(p.length.cat)
  print(p.length.exon)
  # if (length(FL_multisample_indices)>0) {
  #   print(p.length.all.sample)
  #   print(p.length.exon.sample)
  # }
  
  # 2. general parameters by structual categories
  s <- textGrob("Structural Isoform Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p1)
  if (length(p1.s.list) > 0) {
    for (i in 1:length(p1.s.list)) {
      print(p1.s.list[i])
    }
  }
  print(p4)
  print(p4.s1)
  print(p4.s2)
  print(p4.s3)
  print(p5)
  print(p5.s1)
  print(p5.s2)
  print(p5.s3)
  print(pSTM)
  print(pSTM_perc)
  print(pSTM.s1)
  print(pSTM_perc.s1)
  print(pSTM.s2)
  print(pSTM_perc.s2)
  
  if (!all(is.na(data.class$iso_exp))){
    print(p8)
    if (!all(is.na(data.FSMISM$iso_exp))){
      print(p8.s1)
    }
    if (!all(is.na(data.NICNNC$iso_exp))){
      print(p8.s2)
    }
    if (!all(is.na(data.other$iso_exp))){
      print(p8.s3)
    }
  }
  
  
  if (nrow(data.FSM) > 0 ) {
    print(p2)
    print(p3)
  }
  if (!all(is.na(data.class$gene_exp))){
    if (nrow(data.class[data.class$structural_category=="NNC",])!=0){
      print(p12)
    }
  }
  #if (!all(is.na(data.class$gene_exp))){
  #    if (nrow(data.class[data.class$structural_category=="NNC",])!=0 & nrow(data.class[data.class$structural_category=="FSM",])!=0 ){]
  #        print(p13)
  #        print(p13.c)
  #    }
  #}
  
  #3. splice junction
  
  s <- textGrob("Splice Junction Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p23.a)
  print(p23.b)
  #   print(p24)
  #   print(p25)
  #   print(p26)
  #
  #   if (!all(is.na(data.junction$total_coverage)) & !all(is.na(data.class$iso_exp))){
  #     print(pn1.2)
  #   }
  #
  
  if (!all(is.na(data.junction$total_coverage))) {
    print(pn4.a)
    print(pn4.b)
  }
  
  s <- textGrob("Features of Bad Quality", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p28.RTS)
  print(p28.SJ)
  if (n_t3.SJ>0 & n_t3.RTS>0 & !all(is.na(data.class$min_cov))) {
    print(p28.Cov)
  }
  if (n_t3.SJ>0 & n_t3.RTS>0 &!all(is.na(data.class$predicted_NMD))) {
    print(p28.NMD)
  }
  if (n_t3.SJ>0 & n_t3.RTS>0) {
    print(p28)
}
  
  s <- textGrob("Features of Good Quality", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  
  if (!all(is.na(data.class$min_cov))) {
    print(p28.a.Cov)
  }
  
  if (nrow(t3.SJ) > 0) {
    print(p28.a)
  }

  dev.off()
  
  print("SQANTI3 report successfully generated!")
}
  


  
generatePDFreport()

