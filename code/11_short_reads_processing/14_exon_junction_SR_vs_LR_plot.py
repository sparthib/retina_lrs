#!/usr/bin/env python3
# 14_exon_junction_SR_vs_LR_plot.py
# Combine per-sample exon-exon junction distributions (short-read STAR, MAPQ=255)
# with the long-read organoid distributions (code/02_bam_QC output) and plot the
# per-read junction-count distribution by platform. Short reads are physically
# capped at ~2-3 junctions per 150bp read; long reads span many.
import pandas as pd, numpy as np, glob, os
import matplotlib.pyplot as plt

SR_DIR = "/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/star_bams/junction_per_read"
LR_DIR = "/users/sparthib/retina_lrs/processed_data/exon_exon"
LR_RO = ["EP1-BRN3B-RO","EP1-WT_hRO_2","EP1-WT_ROs_D45","H9-BRN3B_hRO_2",
         "H9-BRN3B-RO","H9-CRX_hRO_2","H9-CRX_ROs_D45"]

def load_dist(path):
    return pd.read_csv(path).iloc[0,1:].values.astype(float)

rows=[]
for f in sorted(glob.glob(f"{SR_DIR}/*_junction_per_read.csv")):
    s=os.path.basename(f).replace("_junction_per_read.csv","")
    rows.append(["Short read",s]+list(load_dist(f)))
for s in LR_RO:
    f=f"{LR_DIR}/{s}_junction_per_read.csv"
    if os.path.exists(f):
        rows.append(["Long read",s]+list(load_dist(f)))
wide=pd.DataFrame(rows, columns=["platform","sample"]+[f"junctions_{i}" for i in range(11)])
wide.to_csv("exon_junction_per_read_SR_vs_LR.csv", index=False)

long=wide.melt(id_vars=["platform","sample"], var_name="jc", value_name="frac")
long["junctions"]=long["jc"].str.replace("junctions_","").astype(int)
sub=long[(long.junctions>=1)&(long.junctions<=10)]
sr_col,lr_col="#C0392B","#2166AC"; width=0.38; junc=list(range(1,11))
fig,ax=plt.subplots(figsize=(7.2,4.2))
for k,(plat,col) in enumerate([("Short read",sr_col),("Long read",lr_col)]):
    off=-width/2 if k==0 else width/2
    data=[sub[(sub.platform==plat)&(sub.junctions==j)]["frac"].values*100 for j in junc]
    pos=[j+off for j in junc]
    ax.boxplot(data,positions=pos,widths=width*0.9,patch_artist=True,showfliers=False,
        medianprops=dict(color="black"),boxprops=dict(facecolor=col,alpha=0.25,edgecolor=col),
        whiskerprops=dict(color=col),capprops=dict(color=col))
    for j,d,p in zip(junc,data,pos):
        ax.scatter(np.random.RandomState(j).normal(p,0.04,len(d)),d,s=7,color=col,alpha=0.7,zorder=3)
ax.set_xticks(junc); ax.set_xlabel("Number of exon-exon junctions spanned per read")
ax.set_ylabel("Reads (% of reads in multi-exon PC genes)")
ax.set_title("Long reads span many more exon-exon junctions per read than short reads")
ax.spines[['top','right']].set_visible(False)
fig.tight_layout(); fig.savefig("exon_junction_SR_vs_LR.png", dpi=200)
