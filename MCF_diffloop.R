## loading packages
library(diffloop)
library(diffloopdata)
library(ggplot2)
library(GenomicRanges)
library(ggrepel)
library(DESeq2)
library(dplyr)

## loading FitHiChip Q-value bed files
data1 <- read.delim("FitHiChIP.interactions_MCF7-A_FitHiC_Q0.01.bed")
data2 <- read.delim("FitHiChIP.interactions_MCF7-B_FitHiC_Q0.01.bed")
data3 <- read.delim("FitHiChIP.interactions_MCF10A-A_FitHiC_Q0.01.bed")
data4 <- read.delim("FitHiChIP.interactions_MCF10A-B_FitHiC_Q0.01.bed")

## subsetting the dataframe to required columns
data1.1 <- data1[,c(1:7)]
data2.2 <- data2[,c(1:7)]
data3.3 <- data3[,c(1:7)]
data4.4 <- data4[,c(1:7)]

## Removing the chr character infront of chromosome number
data1.1$chr1 <- gsub("chr","",as.character(data1.1$chr1))
data1.1$chr2 <- gsub("chr","",as.character(data1.1$chr2))
data2.2$chr1 <- gsub("chr","",as.character(data2.2$chr1))
data2.2$chr2 <- gsub("chr","",as.character(data2.2$chr2))
data3.3$chr1 <- gsub("chr","",as.character(data3.3$chr1))
data3.3$chr2 <- gsub("chr","",as.character(data3.3$chr2))
data4.4$chr1 <- gsub("chr","",as.character(data4.4$chr1))
data4.4$chr2 <- gsub("chr","",as.character(data4.4$chr2))

## Inserting dot as column in dataframe
data1.1[,'dot'] <- "."
data2.2[,'dot'] <- "."
data3.3[,'dot'] <- "."
data4.4[,'dot'] <- "."

## Rearranging the dataframe columns
mcf7_1 <- data1.1[,c(1:6,8,7)]
mcf7_2 <- data2.2[,c(1:6,8,7)]
mcf10a_1 <- data3.3[,c(1:6,8,7)]
mcf10a_2 <- data4.4[,c(1:6,8,7)]

mcf7.rep1 <- mcf7_1[(mcf7_1$chr1==15 | mcf7_1$chr2==15),]
mcf7.rep2 <- mcf7_2[(mcf7_2$chr1==15 | mcf7_2$chr2==15),]
mcf10a.rep1 <- mcf10a_1[(mcf10a_1$chr1==15 | mcf10a_1$chr2==15),]
mcf10a.rep2 <- mcf10a_2[(mcf10a_2$chr1==15 | mcf10a_2$chr2==15),]

## Writing bedpe files
write.table(mcf7_1, file = "mcf7_1.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(mcf7_2, file = "mcf7_2.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(mcf10a_1, file = "mcf10a_1.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(mcf10a_2, file = "mcf10a_2.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)


## diffloop downstream analysis with bedpe files
bed_dir <- file.path("/Users/ayalurik/Documents/HiChIP_diffloop/diffloop/")
bed_dir

samples <- c("mcf7_1","mcf7_2","mcf10a_1","mcf10a_2")
full <- loopsMake(bed_dir, samples)
celltypes <- c("mcf7", "mcf7", "mcf10a", "mcf10a")
full <- updateLDGroups(full, celltypes)
head(full, 4)
dim(full)


## Quality control
# remove and loops that merged together from import
full1 <- subsetLoops(full, full@rowData$loopWidth >= 5000)
dim(full1)

# Filtering loops that are FDR > 0.01
noCNV_corrected <- mangoCorrection(full1, FDR = 0.01)
dim(noCNV_corrected)

# filtering loops that are present strongly (>= 5 PETs) in one replicate
# but absent (== 0 PETs) in the other replicate
cm <- noCNV_corrected@counts
k_dis <- ((cm[,1]>=5 &cm[,2]==0)|(cm[,2]>=5&cm[,1]==0))
m_dis <- ((cm[,3]>=5 &cm[,4]==0)|(cm[,4]>=5&cm[,3]==0))
qc_filt <- subsetLoops(noCNV_corrected, !(k_dis | m_dis))
dim(qc_filt)

p1 <- loopDistancePlot(qc_filt)
p1

## counts the number of loops for each sample and returns whether they are single,
## self, unique or none
loopMetrics(qc_filt)

# PC plot
pcp1dat <- qc_filt
pcp1dat@colData$sizeFactor <- 1
samples <- c("MCF7_1", "MCF7_2", "MCF10A_1", "MCF10A_2")
pcp1 <- pcaPlot(pcp1dat) + geom_text_repel(aes(label=samples)) + 
ggtitle("PC Plot with no
  Size Factor Correction") + theme(legend.position="none")
pcp1

samples <- c("MCF7_1", "MCF7_2", "MCF10A_1", "MCF10A_2")
pcp2 <- pcaPlot(qc_filt) + geom_text_repel(aes(label=samples)) + 
  ggtitle("PC Plot with Size Factor Correction") +
  theme(legend.position = "none")
pcp2

## Differential Loop Calling
km_filt <- qc_filt
dim(km_filt)

# First model using edgeR over-dispersed Poisson regression
km_res <- quickAssoc(km_filt)
head(km_res@rowData)

## Epigenetic Annotation
h3dir <- file.path("/Users/ayalurik/Documents/HiChIP_diffloop")
kh3 <- paste0(h3dir, "/", "MCF7_H3k27ac_hg38_narrowpeak.bed")
mh3 <- paste0(h3dir, "/", "MCF10A_H3k27ac_hg38_narrowpeak.bed")

h3k27ac.k <- rmchr(padGRanges(bedToGRanges(kh3), pad = 1000))
h3k27ac.m <- rmchr(padGRanges(bedToGRanges(mh3), pad = 1000))
enhancer <- GenomicRanges::union(h3k27ac.m, h3k27ac.k)
promoter <- padGRanges(getHumanTSS(), pad = 1000)
km_anno <- annotateLoops(km_res, enhancer = enhancer, promoter = promoter)

# writing out results
summary <- as.data.frame(km_anno@anchors)
stats <- as.data.frame(km_anno@rowData)
interactions <- as.data.frame(km_anno@interactions)
data.interactions <- summary(km_anno)
sig_results <- data.interactions[data.interactions$PValue < 0.05 & 
                                   data.interactions$FDR < 0.01, ]
mcf7_up <- data.interactions[data.interactions$logFC >= 1 & data.interactions$FDR <= 0.05, ]
mcf10a_up <- data.interactions[data.interactions$logFC <= -1 & data.interactions$FDR <= 0.05, ]
en_pr <- sig_results[sig_results$loop.type == 'e-p',]
en_en <- sig_results[sig_results$loop.type == 'e-e',]
## saving a rda file
# save(km_anno, file = "diffloop_MCF.rda")

## Visualization
plotTopLoops(km_anno, n = 30, PValue = 0.05, FDR = 0.05, organism = "h",
                   colorLoops = TRUE)
chr1reg <- GRanges(seqnames=c("8"),ranges=IRanges(start=c(128125000),end=c(128450000)))
p1 <- loopPlot(km_anno, chr1reg, organism = 'h')

## multiple loop plots
gr1 <- c(GRanges(seqnames=c("8"),ranges=IRanges(start=c(128125000),end=c(128450000))), 
        GRanges(seqnames=c("1"),ranges=IRanges(start=c(32350000),end=c(32465000))),
        GRanges(seqnames=c("5"),ranges=IRanges(start=c(150370000),end=c(150500000))),
        GRanges(seqnames=c("21"),ranges=IRanges(start=c(41455000),end=c(41815000))),
        GRanges(seqnames=c("3"),ranges=IRanges(start=c(64080000),end=c(64550000))))

manyLoopPlots(km_anno, gr1, organism = 'h', colorLoops = TRUE)

# write.csv(data.interactions, file = "MCF7_vs_MCF10A_diff_interactions.csv")
