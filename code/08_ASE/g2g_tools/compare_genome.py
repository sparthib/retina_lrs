def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a dictionary of sequences with headers.

    Args:
    fasta_file (str): Path to the input FASTA file.

    Returns:
    dict: A dictionary with headers as keys and sequences as values.
    """
    sequences = {}
    header = None
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]  # Remove '>'
                sequences[header] = ""
            elif header:
                sequences[header] += line
    return sequences

def compare_L_R_sequences(fasta_file, output_file):
    """
    Compare sequences labeled with 'L' and 'R' in a genome FASTA file,
    output mismatched sequences, and print a summary of differences.

    Args:
    fasta_file (str): Path to the input FASTA file.
    output_file (str): Path to the output FASTA file with mismatched sequences.
    """
    sequences = parse_fasta(fasta_file)
    lr_pairs = {}
    
    # Organize sequences by chromosome
    for header, sequence in sequences.items():
        if "_L" in header:
            chrom = header.split("_L")[0]
            lr_pairs.setdefault(chrom, {})['L'] = (header, sequence)
        elif "_R" in header:
            chrom = header.split("_R")[0]
            lr_pairs.setdefault(chrom, {})['R'] = (header, sequence)

    # Open the output file to write mismatched sequences
    with open(output_file, 'w') as outfile:
        total = 0
        mismatches = 0
        
        for chrom, lr_seqs in lr_pairs.items():
            if 'L' in lr_seqs and 'R' in lr_seqs:
                total += 1
                l_header, l_seq = lr_seqs['L']
                r_header, r_seq = lr_seqs['R']
                
                # Compare the L and R sequences
                if l_seq != r_seq:
                    mismatches += 1
                    # Write mismatched L and R sequences to the output file
                    outfile.write(f">{l_header}\n{l_seq}\n")
                    outfile.write(f">{r_header}\n{r_seq}\n")

        # Print summary
        print(f"Total L-R pairs compared: {total}")
        print(f"Total mismatched L-R pairs: {mismatches}")
        if total > 0:
            print(f"Proportion of mismatched pairs: {mismatches / total:.2%}")
        else:
            print("No L-R pairs found in the input file.")

if __name__ == "__main__":
    # Path to your genome FASTA file
    fasta_file = "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/g2gtools/SRR1091088/SRR1091088_diploid_genome.fa"
    
    # Output file for sequences with mismatches
    output_file = "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/g2gtools/SRR1091088/mismatched_genome_L_R_sequences.fasta"
    
    compare_L_R_sequences(fasta_file, output_file)
