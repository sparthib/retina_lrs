from Bio import SeqIO

def compare_L_R_sequences(fasta_file, output_file):
    """
    Compare sequences labeled with 'L' and 'R' in a genome FASTA file,
    output mismatched sequences, and print a summary of differences.

    Args:
    fasta_file (str): Path to the input FASTA file.
    output_file (str): Path to the output FASTA file with mismatched sequences.
    """
    sequences = {}
    
    # Read the sequences and store them by their headers
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.id
        seq = str(record.seq)
        
        # Classify as 'L' or 'R' based on the header
        if "_L" in header:
            chrom = header.split("_L")[0]
            sequences.setdefault(chrom, {})['L'] = seq
        elif "_R" in header:
            chrom = header.split("_R")[0]
            sequences.setdefault(chrom, {})['R'] = seq

    # Open the output file to write mismatched sequences
    with open(output_file, 'w') as outfile:
        total = 0
        mismatches = 0
        
        for chrom, lr_seqs in sequences.items():
            if 'L' in lr_seqs and 'R' in lr_seqs:
                total += 1
                # Compare the L and R sequences
                if lr_seqs['L'] != lr_seqs['R']:
                    mismatches += 1
                    # Write mismatched L and R sequences to the output file
                    outfile.write(f">{chrom}_L\n{lr_seqs['L']}\n")
                    outfile.write(f">{chrom}_R\n{lr_seqs['R']}\n")

        # Print summary
        print(f"Total L-R pairs compared: {total}")
        print(f"Total mismatched L-R pairs: {mismatches}")
        if total > 0:
            print(f"Proportion of mismatched pairs: {mismatches / total:.2%}")

if __name__ == "__main__":
    # Path to your genome FASTA file
    fasta_file = "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/g2gtools/SRR1091088/SRR1091088_diploid_genome.fa"
    
    # Output file for sequences with mismatches
    output_file = "mismatched_genome_L_R_sequences.fasta"
    
    compare_L_R_sequences(fasta_file, output_file)
