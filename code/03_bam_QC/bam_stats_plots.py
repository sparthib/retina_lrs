import sys
import os
import pandas as pd
import aplanat
from aplanat.hist import histogram
from bokeh.io import export_png

# Sample, input, and output directories from command line arguments
sample = sys.argv[1]
input_dir = sys.argv[2]

# Read the stats file
df = pd.read_csv(os.path.join(input_dir, sample + ".stats"), sep="\t")

# Create histograms
p1 = histogram(
    [df['read_length']], title="Read lengths",
    x_axis_label="read length / bases", y_axis_label="count")
p1.xaxis.formatter.use_scientific = False

p2 = histogram(
    [df['acc']], title="Read accuracy",
    x_axis_label="% accuracy", y_axis_label="count")

p3 = histogram(
    [df['mapq']], title="MAPQ Score",
    x_axis_label="minimap2 MAPQ score", y_axis_label="count")

# Save each plot separately as PNG
output_dir = sys.argv[3]

# Save read length histogram
filename_p1 = os.path.join(output_dir, sample + "_read_lengths.png")
export_png(p1, filename=filename_p1)

# Save read accuracy histogram
filename_p2 = os.path.join(output_dir, sample + "_read_accuracy.png")
export_png(p2, filename=filename_p2)

# Save MAPQ score histogram
filename_p3 = os.path.join(output_dir, sample + "_mapq_score.png")
export_png(p3, filename=filename_p3)
