#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=read_count
#SBATCH -c 5
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/processed_fastq_read_count.txt
#SBATCH -e logs/processed_fastq_read_count.txt
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

# Output file to store the results
output_file="/dcs04/hicks/data/sparthib/retina_lrs/01_input_fastqs/read_counts.txt"

# Clear the output file if it already exists
> "$output_file"

# Loop through all gzipped fastq files in the specified folder
for file in /dcs04/hicks/data/sparthib/retina_lrs/01_input_fastqs/*.fastq.gz; do
    # Count the number of reads (4 lines per read)
    read_count=$(zcat "$file" | wc -l)
    total_reads=$((read_count / 4))
    
    # Print the filename and total read count to the output file
    echo "$file: $total_reads" >> "$output_file"
done

echo "Read counts saved to $output_file"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"

