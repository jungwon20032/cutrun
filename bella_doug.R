library(DiffBind)
library(tidyverse)
library(edgeR)
library(profileplyr)

# where all the files are
setwd('C:/Users/jungw/thesis')

# name of sample sheet
bellasamples <- read.csv('bellaDoug.csv')

bellasamples
# loading the samples into a dba object
dbObjBella <- dba(sampleSheet=bellasamples)
dbObjBella

# dba object only contains consensus peaks and comparing overlapping peaks
dba.plotHeatmap(dbObjBella)
dba.plotVenn(dbObjBella, dbObjBella$masks$ovary & dbObjBella$masks$High)
dba.plotVenn(dbObjBella, dbObjBella$masks$embryo)


# binding affinity calculations by importing bam files into dba object
dbObjBellaC <- dba.count(dbObjBella, bUseSummarizeOverlaps=TRUE)
dbObjBellaC


# preliminary graphs showing binding affinity stuff
dba.plotHeatmap(dbObjBellaC)
dba.plotPCA(dbObjBellaC)

# this one takes a while and is the heatmap showing region around each peak
profilesBella <- dba.plotProfile(dbObjBellaC)
dba.plotProfile(profilesBella)

# set default analysis threshold for FDA to 0.05
dbObjBellaC$config$th <- 0.05

# for specifying the categories of comparison
categoriesBella <- dba.contrast(dbObjBellaC, minMembers = 2)
categoriesBella


dbObj2 <- dba.peakset(dbObjBella, consensus = -DBA_REPLICATE)
dbObj2 <-  dba(dbObj2, mask=dbObj2$masks$Consensus)
dbObj2
?dba
dba.plotMA(analysisBella, contrast = 3, th = 0.1)
# for doing statistical analysis using DESeq2 and edgeR 
analysisBella <- dba.analyze(categoriesBella)
analysisBella


# some plots to look at the statistical analysis of differentially expressed regions
dba.plotHeatmap(analysisBella, contrast = 3, th = 0.1)
dba.plotVenn(analysisBella, contrast=1, method = DBA_ALL_METHODS)
dba.plotMA(analysisBella, contrast = 3, bFlip = TRUE)
dba.plotVolcano(analysisBella, bFlip = TRUE, contrast = 3, th = 0.1)
dba.plotBox(analysisBella)

# for putting the analysis into a GRanges object (th variable to adjust threshold)
reportBella <- dba.report(analysisBella, contrast = 1, bCounts = TRUE, th = 0.05)
reportBella
hist(reportBella$"p-value")

# to transform into a data frame that is more manageable
outputBella <- as.data.frame(reportBella)


# to filter significant peaks that are on x chromosome
xchrom <- outputBella %>%
  filter(seqnames == 'X')

# to transform enriched areas into a bed file that can be visualized on IGV
embryo_enrichac <- outputBella %>%
  filter(FDR < 0.05 & Fold > 0) %>%
  select(seqnames, start, end)
write.table(embryo_enrichac, file = 'embryo_enriched_ac.bed', sep='\t', quote=F, row.names=F, col.names=F)

EAD_enrichac <- outputBella %>%
  filter(FDR < 0.05 & Fold < 0) %>%
  select(seqnames, start, end)
write.table(EAD_enrichac, file = 'EAD_enriched_ac.bed', sep='\t', quote=F, row.names=F, col.names=F)

# Biomart stuff and ChIP Seeker###########################
library(biomaRt)
library(ChIPpeakAnno)

ensembl <- useEnsembl(biomart = 'genes')
searchDatasets(mart = ensembl, pattern = 'melanogaster')

ensembl <- useDataset(dataset = 'dmelanogaster_gene_ensembl', mart = ensembl)
TSS <- getAnnotation(
  ensembl,
  featureType = 'TSS'
)

rep <- dba.report(analysisBella, contrast = 2, th = 0.1)
repoutput <- as.data.frame(rep)
data.peaksAnno = annotatePeakInBatch(rep, AnnotationData = TSS)

head(data.peaksAnno)
gene_annotations <- as.data.frame(data.peaksAnno)

gene_names <- gene_annotations %>%
  filter(startsWith(feature, "FBgn"))

write.csv(gene_names, "bellagenes_High_Zero_DOUG0.1.csv")
