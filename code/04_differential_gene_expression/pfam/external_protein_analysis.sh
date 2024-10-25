#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 20
#SBATCH --job-name=protein_analysis
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/RO_protein_analysis.txt
#SBATCH -e logs/RO_protein_analysis.txt
#SBATCH --time=2-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

# comparisons=("FT_vs_RGC" "RO_D100_vs_RO_D45" "RO_D100_vs_RO_D200" "RO_D200_vs_RO_D45")
comparisons=("ROs")

INPUT_DIR=/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/$item/fastas
OUTPUT_DIR=/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/$item/external_protein_analyses

source activate CPC2
cd $CPC_HOME
for item in ${comparisons[@]}; do
    echo $item
    NT_FASTA=$INPUT_DIR/isoformSwitchAnalyzeR_isoform_nt.fasta
    CPC2_OUTPUT=$OUTPUT_DIR/CPC2_output
    bin/CPC2.py -i $NT_FASTA -o $CPC2_OUTPUT
done

conda deactivate


#### pfam
# https://vcru.wisc.edu/simonlab/bioinformatics/programs/install/pfamscan.htm

source activate pfam
cd /dcs04/hicks/data/sparthib/retina_lrs/PfamScan/pfam_scan
for item in ${comparisons[@]}; do
    echo $item
    ./pfam_scan.py $INPUT_DIR/isoformSwitchAnalyzeR_isoform_AA.fasta \
    ../ -out $OUTPUT_DIR/pfam_results.csv -cpu $SLURM_CPUS_PER_TASK
done
conda deactivate

#### SignalP ####

source activate SignalP
for item in ${comparisons[@]}; do
    echo $item
    signalp6 --fastafile $INPUT_DIR/isoformSwitchAnalyzeR_isoform_AA.fasta \
    --organism eukarya --output_dir $OUTPUT_DIR --format txt --mode fast
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
