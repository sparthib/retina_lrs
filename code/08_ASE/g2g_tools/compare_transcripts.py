# Script to compare transcript sequences, calculate proportions of differences, and output to a new FASTA file

def compare_sequences(seq1, seq2):
    """
    Compare two sequences and identify base pair differences.

    Args:
    seq1 (str): The first sequence (e.g., transcript L).
    seq2 (str): The second sequence (e.g., transcript R).

    Returns:
    list: A list of tuples with the position and the differing base pairs.
    """
    differences = []
    length = min(len(seq1), len(seq2))  # Compare only up to the length of the shorter sequence
    
    for i in range(length):
        if seq1[i] != seq2[i]:
            differences.append((i + 1, seq1[i], seq2[i]))  # Append 1-based index
    
    # Handle the case where one sequence is longer than the other
    if len(seq1) != len(seq2):
        longer_seq = seq1 if len(seq1) > len(seq2) else seq2
        for i in range(length, len(longer_seq)):
            differences.append((i + 1, longer_seq[i], '-'))  # Mark the extra bases in the longer sequence

    return differences

def read_fasta(filename):
    """
    Read sequences from a FASTA-like file.

    Args:
    filename (str): The path to the FASTA file.

    Returns:
    dict: A dictionary where keys are the transcript names and values are the sequences.
    """
    sequences = {}
    with open(filename, 'r') as file:
        current_seq_name = None
        current_seq = []
        
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_seq_name:
                    sequences[current_seq_name] = ''.join(current_seq)
                current_seq_name = line[1:]  # Save transcript name without ">"
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_seq_name:
            sequences[current_seq_name] = ''.join(current_seq)
    
    return sequences

def write_fasta(output_file, name, sequence):
    """
    Write a sequence to a FASTA file.

    Args:
    output_file (file object): Opened file object to write the sequence.
    name (str): The name of the sequence.
    sequence (str): The sequence to write.
    """
    output_file.write(f">{name}\n")
    output_file.write(f"{sequence}\n")

def find_differences(transcripts_file, output_fasta):
    """
    Compare all R and L transcripts from the file, print the proportion of differences, 
    and write the differing sequences to a new FASTA file.

    Args:
    transcripts_file (str): Path to the file containing the transcript sequences.
    output_fasta (str): Path to the output FASTA file for sequences with differences.
    """
    sequences = read_fasta(transcripts_file)
    total_pairs = 0
    differing_pairs = 0
    
    with open(output_fasta, 'w') as out_fasta:
        for name, seq in sequences.items():
            if name.endswith('_L'):
                # Get the corresponding R transcript name
                transcript_R_name = name[:-2] + '_R'
                if transcript_R_name in sequences:
                    seq_R = sequences[transcript_R_name]
                    total_pairs += 1
                    differences = compare_sequences(seq, seq_R)
                    
                    if differences:
                        differing_pairs += 1
                        print(f"Differences between {name} and {transcript_R_name}:")
                        for pos, base_L, base_R in differences:
                            print(f"Position {pos}: {base_L} (L) vs {base_R} (R)")
                        
                        # Write both sequences to the output FASTA file
                        write_fasta(out_fasta, name, seq)
                        write_fasta(out_fasta, transcript_R_name, seq_R)
                else:
                    print(f"Transcript {transcript_R_name} not found.")
        
        # Calculate and print proportion of differing pairs
        if total_pairs > 0:
            proportion_differing = differing_pairs / total_pairs
            print(f"\nTotal transcript pairs: {total_pairs}")
            print(f"Transcript pairs with differences: {differing_pairs}")
            print(f"Proportion of pairs with differences: {proportion_differing:.2f}")
        else:
            print("No L and R transcript pairs found.")

if __name__ == "__main__":
    # Replace with the path to your transcript file
    transcripts_file = "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/g2gtools/SRR1091088/extracted_transcripts.fasta"
    
    # Replace with the path to your output FASTA file
    output_fasta = "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/g2gtools/SRR1091088/differing_transcripts.fasta"
    
    find_differences(transcripts_file, output_fasta)
    

    
    
    
    
