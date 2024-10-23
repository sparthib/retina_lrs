import re

def extract_transcripts(log_file, fasta_output):
    """
    Extract transcript names and sequences from a log file and save them to a FASTA file.

    Args:
    log_file (str): Path to the log file containing the transcript information.
    fasta_output (str): Path to the output FASTA file.
    """
    with open(log_file, 'r') as infile, open(fasta_output, 'w') as outfile:
        current_transcript = None
        current_sequence = []
        
        for line in infile:
            line = line.strip()

            # Match transcript names in the format >ENST00000..._L or >ENST00000..._R
            if re.match(r'>ENST\d+.\d+_(L|R)', line):
                if current_transcript:
                    # Write the previous transcript and sequence to the FASTA file
                    outfile.write(f">{current_transcript}\n")
                    outfile.write(f"{''.join(current_sequence)}\n")
                
                # Start a new transcript
                current_transcript = line[1:]  # Remove the ">" from the name
                current_sequence = []
            elif not line.startswith('*') and not re.match(r'\d{4}-\d{2}-\d{2}', line):
                # Add sequence lines (ignore job info and separators like "****")
                current_sequence.append(line)

        # Write the last transcript and sequence
        if current_transcript:
            outfile.write(f">{current_transcript}\n")
            outfile.write(f"{''.join(current_sequence)}\n")

if __name__ == "__main__":
    # Path to your log file
    log_file = "/users/sparthib/retina_lrs/code/08_ASE/g2g_tools/logs/g2gtools.txt"
    
    # Path to the output FASTA file
    fasta_output = "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/g2gtools/SRR1091088/extracted_transcripts.fasta"
    
    extract_transcripts(log_file, fasta_output)
