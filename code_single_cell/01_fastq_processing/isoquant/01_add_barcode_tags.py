import pysam
import argparse

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Add CB and UMI tags to BAM reads.")
parser.add_argument("input_bam", help="Input BAM file")
parser.add_argument("output_bam", help="Output BAM file with tags added")
args = parser.parse_args()

# Open the BAM file
input_bam = args.input_bam
output_bam = args.output_bam  

# Open the input BAM file in read mode
bamfile = pysam.AlignmentFile(input_bam, "rb")
# Open a new BAM file for writing, using the same header as the input
output_bamfile = pysam.AlignmentFile(output_bam, "wb", header=bamfile.header)

for read in bamfile:
    # Get the read name (e.g., CATTGAGGTGATGAAT_AAATCTGCTTGT#109dbc01-6a19-4ee7-9d45-1a42a1584da7_-)
    read_name = read.query_name
    
    # Split the read name into barcode and UMI
    barcode, rest = read_name.split('_', 1)
    umi = rest.split('#')[0]  # UMI is before the '#'

    # Add the CB (cell barcode) and UMI tags
    read.set_tag("CB", barcode, value_type='Z')  # 'Z' indicates a string tag
    read.set_tag("UM", umi, value_type='Z')

    # Write the modified read to the new BAM file
    output_bamfile.write(read)

# Close the files
bamfile.close()
output_bamfile.close()
