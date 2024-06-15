import argparse
import gzip
from collections import defaultdict
from Bio import SeqIO

def read_fasta(file_path):
    fasta_dict = defaultdict(list)
    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_dict[record.id].append(str(record.seq))
    return fasta_dict

def compare_sequences(fasta_dict):
    discrepancies = {}
    duplicate_counts = []

    for read_id, sequences in fasta_dict.items():
        unique_sequences = set(sequences)
        if len(unique_sequences) > 1:
            discrepancies[read_id] = sequences
        if len(sequences) > 1:
            duplicate_counts.append(len(sequences))
    
    total_reads = len(fasta_dict)
    total_duplicates = len(duplicate_counts)
    proportion_duplicates = total_duplicates / total_reads if total_reads > 0 else 0

    return discrepancies, duplicate_counts, total_reads, total_duplicates, proportion_duplicates

def main():
    parser = argparse.ArgumentParser(description="Compare FASTA sequences with the same read ID.")
    parser.add_argument("file_path", help="Path to the gzipped FASTA file.")
    args = parser.parse_args()

    fasta_dict = read_fasta(args.file_path)
    discrepancies, duplicate_counts, total_reads, total_duplicates, proportion_duplicates = compare_sequences(fasta_dict)
    
    if discrepancies:
        print("Discrepancies found:")
        for read_id, sequences in discrepancies.items():
            print(f"Read ID: {read_id}")
            for seq in sequences:
                print(f"  Sequence: {seq}")
    else:
        print("All sequences with the same read ID are identical.")

    print("\nSummary Statistics:")
    print(f"Total reads: {total_reads}")
    print(f"Total duplicates: {total_duplicates}")
    print(f"Proportion of reads that are duplicates: {proportion_duplicates:.2%}")

    if duplicate_counts:
        print(f"Maximum duplicates for a single read ID: {max(duplicate_counts)}")
        print(f"Average duplicates per read ID: {sum(duplicate_counts) / len(duplicate_counts):.2f}")
    else:
        print("No duplicate reads found.")

if __name__ == "__main__":
    main()

