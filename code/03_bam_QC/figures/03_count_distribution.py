import argparse
from collections import Counter

# Function to process the input file and extract counts
def process_count_file(input_file):
    counts = []
    
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()  # Remove extra whitespace/newlines
            # Split the line by comma and extract the count (second element)
            try:
                _, count = line.split(',')
                count = int(count)  # Convert the count to an integer
                counts.append(count)
            except ValueError:
                print(f"Invalid line or format: {line}")
    
    # Count the frequency of each unique count
    count_frequencies = Counter(counts)
    
    return count_frequencies

# Function to write the unique counts and their frequencies to an output file
def write_count_frequencies_to_file(output_file, count_frequencies):
    with open(output_file, 'w') as out_f:
        out_f.write("Unique_Count,Frequency\n")
        for count, freq in sorted(count_frequencies.items()):
            out_f.write(f"{count},{freq}\n")

    print(f"Count frequencies written to {output_file}")

# Main function to handle argument parsing
def main():
    parser = argparse.ArgumentParser(description="Compute frequency distribution of counts from input file")
    
    # Define the input argument for the sample file
    parser.add_argument("input_file", help="Path to the input file containing CIGAR strings and counts")
    
    # Optional argument for the output file (default name if not provided)
    parser.add_argument("-o", "--output_file", default="count_distribution_output.txt",
                        help="Path to the output text file (default: count_distribution_output.txt)")
    
    args = parser.parse_args()
    
    # Process the file to get count frequencies
    count_frequencies = process_count_file(args.input_file)
    
    # Write the results to the output file
    write_count_frequencies_to_file(args.output_file, count_frequencies)

if __name__ == "__main__":
    main()
