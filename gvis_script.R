library(GenomicFeatures)
library(Gviz)
library(viridis)
library(scales)

# This script imports signal tracks and displays them at the specified regions of the human 
# genome, listed at the end of the script.

cutntag_dir <- file.path("fly_data")
gviz_dir <- file.path("gviz_dir")
# Create genomic axis track
axis.track <- GenomeAxisTrack(col = "black", scale = 0.1, col.range = "black")

options(ucscChromosomeNames = FALSE)

# Create gene annotation track
txdb <- makeTxDbFromGFF(file.path("fly_data", "dros_annot.gtf"))

genome.track <- GeneRegionTrack(txdb, genome = "FG", shape = "arrow", names = "Genes", col = "black",
                                showId = TRUE, fill = "black", trancriptAnnotation = "gene_symbol")

# Create data tracks

# Get colors from the viridis palette for the number of tracks to be plotted
show_col(viridis_pal()(20))

# Set y-axis limits

h3k27accnt_lim <- c(0,50)
h3k27me3cnt_lim <- c(0,0.5)

# K562 histone marks 
H3K27ac <- DataTrack(range = file.path(cutntag_dir, "H3K27ac_merged.spikenorm.bw"), genome = "fly_data",
                     name = "H3K27ac CUT&RUN", col.histogram = "#3CBC75FF", fill.histogram = "#3CBC75FF", ylim = h3k27accnt_lim)

#H3K4me1 <- DataTrack(range = file.path(cutntag_dir, "H3K4me1_merged.spikenorm.bw"), genome = "GRCh38",
 #                    name = "H3K4me1 CUT&Tag", col.histogram = "#39558CFF", fill.histogram = "#39558CFF", ylim = h3k4cnt_lim) 

#H3K4me2 <- DataTrack(range = file.path(cutntag_dir, "H3K4me2_merged.spikenorm.bw"), genome = "GRCh38",
#                     name = "H3K4me2 CUT&Tag", col.histogram = "#287D8EFF", fill.histogram = "#287D8EFF", ylim = h3k4cnt_lim) 

#H3K4me3 <- DataTrack(range = file.path(cutntag_dir, "H3K4me3_merged.spikenorm.bw"), genome = "GRCh38",
#                     name = "H3K4me3 CUT&Tag", col.histogram = "#1F968BFF", fill.histogram = "#1F968BFF", ylim = h3k4cnt_lim)

#IgG <- DataTrack(range = file.path(cutntag_dir, "IgG_merged.spikenorm.bw"), genome = "GRCh38",
#                 name = "IgG CUT&Tag", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = h3k27cnt_lim) 

# K562 enhancer annotations
enhancers <- AnnotationTrack(range = file.path("fly_data", "redfly.bed"), genome = "fly_data",
                             name = "Enhancers", col = "black", col.line = "white", fill = "black")

# STRIPE-seq - TBPL1
cairo_pdf(file = file.path(gviz_dir, "eyeless.pdf"), height = 12, width = 12)
plotTracks(list(axis.track, 
                H3K27ac,
                genome.track), 
           chromosome = "4", from = 697689, to = 721173, 
           background.title = "white", 
           col.title = "black", 
           col.axis = "black", 
           type = "histogram", 
           baseline = 0, 
           col.baseline = "black")
dev.off()
