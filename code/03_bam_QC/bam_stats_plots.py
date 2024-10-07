# use python to plot the results
import sys
import os
import pandas as pd
import aplanat
from aplanat.hist import histogram
from bokeh.layouts import gridplot

sample = sys.argv[1]
df = pd.read_csv(sample + ".alignment.bam.stats", sep="\t")
p1 = histogram(
    [df['read_length']], title="Read lengths",
    x_axis_label="read length / bases", y_axis_label="count")
p1.xaxis.formatter.use_scientific = False

p2 = histogram(
    [df['acc']], title="Read accuracy",
    x_axis_label="% accuracy", y_axis_label="count")

p3=histogram(
    [df['mapq']], title="MAPQ Score",
    x_axis_label="minimap2 MAPQ score ", y_axis_label="count")
    
grid = gridplot((p1, p2, p3), ncols=2)
# save the plot
aplanat.save(grid, filename=sample + "_bam_stats.pdf")

