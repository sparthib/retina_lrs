from Bio import SeqIO

def read_fasta(file_path):
    fasta_dict = {}
    for record in SeqIO.parse(file_path, "fasta"):
        if record.id in fasta_dict:
            fasta_dict[record.id].append(str(record.seq))
        else:
            fasta_dict[record.id] = [str(record.seq)]
    return fasta_dict

def compare_sequences(fasta_dict):
    discrepancies = {}
    for read_id, sequences in fasta_dict.items():
        if len(set(sequences)) > 1:
            discrepancies[read_id] = sequences
    return discrepancies

def main(file_path):
    fasta_dict = read_fasta(file_path)
    discrepancies = compare_sequences(fasta_dict)
    if discrepancies:
        print("Discrepancies found:")
        for read_id, sequences in discrepancies.items():
            print(f"Read ID: {read_id}")
            for seq in sequences:
                print(f"  Sequence: {seq}")
    else:
        print("All sequences with the same read ID are identical.")

if __name__ == "__main__":
    file_path = "your_fasta_file.fasta"
    main(file_path)
