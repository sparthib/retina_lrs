import pysam
import argparse

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Add CB and UMI tags to BAM reads, remove duplicates, and keep the longest read for each UMI.")
parser.add_argument("input_bam", help="Input BAM file")
parser.add_argument("output_bam", help="Output BAM file with tags added")
args = parser.parse_args()

# Open the input and output BAM files
input_bam = args.input_bam
output_bam = args.output_bam

# Open the input BAM file in read mode
bamfile = pysam.AlignmentFile(input_bam, "rb")
# Open a new BAM file for writing, using the same header as the input
output_bamfile = pysam.AlignmentFile(output_bam, "wb", header=bamfile.header)

# Dictionary to store the longest read for each (barcode, UMI) pair
umi_dict = {}

for read in bamfile:
    # Get the read name (e.g., CATTGAGGTGATGAAT_AAATCTGCTTGT#109dbc01-6a19-4ee7-9d45-1a42a1584da7_-)
    read_name = read.query_name
    
    # Split the read name into barcode and UMI
    barcode, rest = read_name.split('_', 1)
    umi = rest.split('#')[0]  # UMI is before the '#'

    # Get the current read length
    read_length = read.infer_read_length()

    # Check if this (barcode, UMI) pair is already in the dictionary
    if (barcode, umi) in umi_dict:
        # Compare lengths, keep the longest read
        if read_length > umi_dict[(barcode, umi)].infer_read_length():
            umi_dict[(barcode, umi)] = read  # Update with the longer read
    else:
        # Add this (barcode, UMI) pair to the dictionary
        umi_dict[(barcode, umi)] = read

# Write the longest reads for each (barcode, UMI) to the output BAM file
for read in umi_dict.values():
    # Add the CB (cell barcode) and UMI tags
    barcode, rest = read.query_name.split('_', 1)
    umi = rest.split('#')[0]
    read.set_tag("CB", barcode, value_type='Z')  # 'Z' indicates a string tag
    read.set_tag("UM", umi, value_type='Z')

    # Write the modified read to the new BAM file
    output_bamfile.write(read)

# Close the files
bamfile.close()
output_bamfile.close()
