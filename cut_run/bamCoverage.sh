#PBS -N coverage
#PBS -l procs=20
#PBS -l mem=20gb
#PBS -q una
#PBS -M jlee3@hamilton.edu
#PBS -m ea
#PBS -j oe
#PBS -l walltime=12:00:00
#PBS -r n

source /home/jlee3/ChIPseq/settings.conf

#########################################
## ChIP-seq Analysis of Multiple Samples
#########################################

## loading conda environment
## -------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq-automation

# Directory containing E. coli aligned bam files
ecoli_aligned_dir="/home/jlee3/ChIPseq/results/ecoli_aligned"

# Directory containing Drosophila aligned bam files
drosophila_aligned_dir="/home/jlee3/ChIPseq/results/aligned"

# Normalize Drosophila aligned bam files and generate bigwig files for each E. coli sample
for ecoli_bam_file in "${ecoli_aligned_dir}"/*.bam; do
    # Calculate scale factor for each E. coli sample
    ecoli_reads=$(samtools view -c -F 4 "${ecoli_bam_file}")
    scale_factor=$(bc <<< "scale=3; 10000 / ${ecoli_reads}")

    # Normalize Drosophila aligned bam files and generate bigwig files
    for drosophila_bam_file in "${drosophila_aligned_dir}"/*.bam; do
        output_bigwig="${drosophila_bam_file%.bam}_normalized.bigwig"
        bamCoverage --bam "${drosophila_bam_file}" -o "${output_bigwig}" --scaleFactor "${scale_factor}"
    done
done
