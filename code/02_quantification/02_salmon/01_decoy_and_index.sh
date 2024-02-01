#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=create_index
#SBATCH -c 2
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/index.txt
#SBATCH -e logs/index.txt


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"


source activate salmon 

decoy_location=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/salmon/decoy
transcriptome_fa=/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts_short_header.fa
genome_fa=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
gentrome_folder=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/salmon/gentrome/
salmon_index=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/salmon/transcripts/index

# mkdir -p $decoy_location
# mkdir -p $gentrome_folder
# 
# grep "^>" $genome_fa | cut -d " " -f 1 > $decoy_location/decoys.txt
# 
# sed -i.bak -e 's/>//g' $decoy_location/decoys.txt
# 
# cat $transcriptome_fa $genome_fa > $gentrome_folder/gentrome.fa


salmon index -t $gentrome_folder/gentrome.fa -d $decoy_location/decoys.txt \
          -i $salmon_index  -p 12  --gencode
          
conda deactivate salmon

echo "**** Job ends ****"
date +"%Y-%m-%d %T"