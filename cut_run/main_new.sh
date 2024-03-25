#PBS -N normalization
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

## bowtie 2 read alignment
## -----------------------
# Create directory for E. coli genomic index
mkdir -p ${OUTDIR}/ecoli_genome/index

# Create Bowtie2 index for E. coli genome
bowtie2-build \
-f --threads $CORES \
/home/jlee3/ChIPseq/normalization/ecoli.fa \
${OUTDIR}/ecoli_genome/index/ecoli_genome

# Create directory to output alignments for E. coli genome
mkdir -p ${OUTDIR}/ecoli_aligned

# Read alignment against E. coli genome
Rscript ${BASEDIR}/bin/ecolialignReads.R \
--outdir ${OUTDIR} \
--seqdir $SEQDIR \
--threads $CORES \
--samplesheet $SAMPLE_SHEET \

# Convert aligned SAM files to coordinate sorted BAM with index for E. coli genome
for SAM in ${OUTDIR}/ecoli_aligned/*.sam; do
    samtools sort \
    -O BAM -@ $CORES \
    -o ${OUTDIR}/ecoli_aligned/$(basename $SAM .sam).bam \
    $SAM
    samtools index ${OUTDIR}/ecoli_aligned/$(basename $SAM .sam).bam
done

# Index BAM files for E. coli genome
for BAM in ${OUTDIR}/ecoli_aligned/*.bam; do
    samtools index $BAM
