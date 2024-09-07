import argparse

# Function to process the CIGAR file and count 'N's
def count_Ns_in_cigar_file(input_file, output_file):
    with open(output_file, 'w') as out_f:
        out_f.write("ReadLine,Number_of_Ns\n")
        
        with open(input_file, 'r') as file:
            for line in file:
                line = line.strip()  # Remove extra whitespace/newlines
                count = line.count('N')  # Count occurrences of 'N'
                out_f.write(f"{line},{count}\n")

    print(f"Output written to {output_file}")
    
# Main function to handle argument parsing
def main():
    parser = argparse.ArgumentParser(description="Count 'N's in CIGAR strings")
    
    # Define the input argument for the sample file
    parser.add_argument("input_file", help="Path to the input CIGAR file")
    
    # Optional argument for the output file (default name if not provided)
    parser.add_argument("-o", "--output_file", default="cigar_N_counts_output.txt",
                        help="Path to the output text file (default: cigar_N_counts_output.txt)")
    
    args = parser.parse_args()
    
    # Call the function to process the file
    count_Ns_in_cigar_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
    
    
