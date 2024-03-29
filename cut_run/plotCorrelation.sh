#PBS -N correlation
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

## making the plotCorrelatiion directory
mkdir -p ${OUTDIR}/plotCorrelation

multiBamSummary bins --bamfiles ${OUTDIR}/aligned/*.bam -o ${OUTDIR}/plotCorrelation/XYZ.npz

plotCorrelation -in ${OUTDIR}/plotCorrelation/XYZ.npz --whatToPlot heatmap --corMethod spearman -o ${OUTDIR}/plotCorrelation/XYZ.png --labels $(ls ${OUTDIR}/aligned/*.bam | xargs -n 1 basename | sed 's/\.bam//g')
