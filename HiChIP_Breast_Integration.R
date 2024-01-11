source('~/Documents/Rscripts/Hichip_helper_functions.R')


## Load MCF10A and MCF7 interaction files ##

mcf10a_1 <- fread("~/Documents/Breast_Hichip/FitHiChIP.interactions_MCF10A-A_FitHiC_Q0.01.bed")
mcf10a_2 <- fread("~/Documents/Breast_Hichip/FitHiChIP.interactions_MCF10A-B_FitHiC_Q0.01.bed")

mcf7_1 <- fread("~/Documents/Breast_Hichip/FitHiChIP.interactions_MCF7-A_FitHiC_Q0.01.bed")
mcf7_2 <- fread("~/Documents/Breast_Hichip/FitHiChIP.interactions_MCF7-B_FitHiC_Q0.01.bed")


#OverlapLoop() from HiChiP_Overlap_reps.R script
mcf10a_ovp <- OverlapLoop(as.data.frame(mcf10a_1), as.data.frame(mcf10a_2), boundary = 1000, offset = 5000)
mcf10a_ovp$A_AND_B.df %>% dim

mcf7_ovp <- OverlapLoop(as.data.frame(mcf7_1), as.data.frame(mcf7_2), boundary = 1000, offset = 5000)
##


## prep data for the functions above ##
mcf10a <-as.data.frame(mcf10a_2[,c("chr1", "s1", "e1", "chr2", "s2", "e2", "cc", "Q-Value_Bias")])
write.table(mcf10a, file = "~/Documents/Breast_Hichip/MCF10A.interactions.all.mango", row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t")
mcf7 <-as.data.frame(mcf7_2[,c("chr1", "s1", "e1", "chr2", "s2", "e2", "cc", "Q-Value_Bias")])
write.table(mcf7, file = "~/Documents/Breast_Hichip/MCF7.interactions.all.mango", row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t")

mcf10a = loopsMake.mango.null('~/Documents/Breast_Hichip', 'MCF10A',
                              mergegap=0, pad=5000, reduce=FALSE)
mcf10a <- rmchr(mcf10a)
mcf7 = loopsMake.mango.null('~/Documents/Breast_Hichip', samples = 'MCF7',
                              mergegap=0, pad=5000, reduce=FALSE)
mcf7 <- rmchr(mcf7)
##                                  ##

###                   Get TSS for hg19              ###
library(EnsDb.Hsapiens.v75)
library(GenomicDistributionsData)
TSS <- GenomicDistributionsData::TSS_hg19()
promoters <- TSS + 500
promoters <- rmchr(promoters)
##                                                  ###

##                H3k27ac peaks             ###
k27.mcf10a <- fread("~/Documents/Breast_Hichip/MCF10A_H3k27ac_hg38_narrowpeak.bed")
k27.mcf10a$h3k27peakid = paste(k27.mcf10a$V1, k27.mcf10a$V2, k27.mcf10a$V3, sep="_")
k27.mcf10a <- makeGRangesFromDataFrame(k27.mcf10a, keep.extra.columns = T,
                                       seqnames.field="V1", start.field="V2", end.field="V3")
k27.mcf10a = rmchr(k27.mcf10a)

k27.mcf7 <- fread("~/Documents/Breast_Hichip/MCF7_H3k27ac_hg38_narrowpeak.bed")
k27.mcf7$h3k27peakid = paste(k27.mcf7$V1, k27.mcf7$V2, k27.mcf7$V3, sep="_")
k27.mcf7 <- makeGRangesFromDataFrame(k27.mcf7, keep.extra.columns = T,
                                       seqnames.field="V1", start.field="V2", end.field="V3")
k27.mcf7 <- rmchr(k27.mcf7)
##                                            ###


##                Enhancers                    ##
enhancers_mcf10a <- k27.mcf10a[!(k27.mcf10a %over% promoters)]
#enhancers_mcf7<- k27.mcf7[!(k27.mcf7 %over% promoters)]

##                                            ##


###                   MCF7 ChromHMM calls from public data              ###
mcf7_hmm <- fread("~/Documents/Breast_Hichip/GSE57498_MCF7_ChromHMM.bed.gz", skip = 1)
mcf7_hmm$Enhancerpeakid <- paste(mcf7_hmm$V1, mcf7_hmm$V2, mcf7_hmm$V3, sep="_")
mcf7_hmm <- makeGRangesFromDataFrame(mcf7_hmm, keep.extra.columns = T,
                                     seqnames.field="V1", start.field="V2", end.field="V3")
mcf7_enhancer <- mcf7_hmm[mcf7_hmm$V4 == "Enhancer"]
mcf7_enhancer <- rmchr(mcf7_enhancer)

##                                                                      ###


##                            FINE MAPPING/CCV LIST PREP FILES                        ##
brca.fm <- read_excel("~/Documents/Breast_Hichip/41588_2019_537_MOESM3_ESM.xlsx", sheet = 1, skip = 2)
brca.ccv <- read_excel("~/Documents/Breast_Hichip/41588_2019_537_MOESM4_ESM (1).xlsx", sheet = 3, skip = 2)
brca.ccv$end <- brca.ccv$position + 1
brca.ccv$ccvid <- paste0(brca.ccv$chr,brca.ccv$position, brca.ccv$end, sep="_" )
brca.ccv.gr <- makeGRangesFromDataFrame(brca.ccv, seqnames.field="chr", start.field="position",
                                        end.field="end", keep.extra.columns = T)

brca.ccv.gr <- rmchr(brca.ccv.gr)
###                                                                                   ##

###                     TWAS                          ###
brca.twas <- fread("~/Documents/Breast_Hichip/BRCA_TWAS.txt", header = F)
brca.twas <- brca.twas %>% separate(col = "V2", into = c("chr", "rest"), sep = ":" )
brca.twas <- brca.twas %>% separate(col = "rest", into = c("start", "end"), sep = "–" )
brca.twas[,2:4] <- lapply(brca.twas[,2:4], as.numeric)
# 
tmp <- brca.twas[which((brca.twas$end - brca.twas$start) <0) , ]$start
brca.twas[which((brca.twas$end - brca.twas$start) <0) , ]$start <- brca.twas[which((brca.twas$end - brca.twas$start) <0) , ]$end
brca.twas[which((brca.twas$end - brca.twas$start) <0) , ]$end <- tmp
#
colnames(brca.twas)[1] <- "Gene_TWAS"
brca.twas <- makeGRangesFromDataFrame(brca.twas, keep.extra.columns = T)
###                                                   ###


###                   eQTL                              ###
brca.eqtl <- fread("~/Documents/Breast_Hichip/BRCA_eQTL_genes.txt", header = F)
# need to get gene locations from symbol
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
bm <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
      filters=c('hgnc_symbol'),
      values=brca.eqtl$V1,
      mart=ensembl)
bm <- bm[-c(grep('PATCH|CTG', bm$chromosome_name)),]
colnames(bm)[4] <- 'eQTL'
brca.eqtl <- makeGRangesFromDataFrame(bm, seqnames.field ="chromosome_name",
                                      start.field="start_position",
                                      end.field="end_position", keep.extra.columns = T)
###                                                     ###


###                     Differential Looping                  ###
## Negative fold change - enriched in MCF10A; positive - enriched in MCF7 - Important!!
diffloop.breast <- fread('~/Documents/Breast_Hichip/DiffLoop_Breast_ResultTable.txt')
diffloop.mcf10a <- diffloop.breast[diffloop.breast$logFC <= -1 & diffloop.breast$FDR <= 0.05, ]
diffloop.mcf7 <- diffloop.breast[diffloop.breast$logFC >= 1 & diffloop.breast$FDR <= 0.05, ]

#make granges
diffloop.mcf10a_1 <- makeGRangesFromDataFrame(diffloop.mcf10a, seqnames.field="chr_1", start.field="start_1",
                                              end.field="end_1", keep.extra.columns = T)
diffloop.mcf10a_1$DiffMCF10.anchor_1 <- paste(seqnames(diffloop.mcf10a_1), start(diffloop.mcf10a_1), 
                                              end(diffloop.mcf10a_1), sep="_")
diffloop.mcf10a_2 <- makeGRangesFromDataFrame(diffloop.mcf10a, seqnames.field="chr_2", start.field="start_2",
                                              end.field="end_2", keep.extra.columns = T)
diffloop.mcf10a_2$DiffMCF10.anchor_2 <- paste(seqnames(diffloop.mcf10a_2), start(diffloop.mcf10a_2), 
                                              end(diffloop.mcf10a_2), sep="_")

diffloop.mcf7_1 <- makeGRangesFromDataFrame(diffloop.mcf7, seqnames.field="chr_1", start.field="start_1",
                                              end.field="end_1", keep.extra.columns = T)
diffloop.mcf7_1$DiffMC7.anchor_1 <- paste(seqnames(diffloop.mcf7_1), start(diffloop.mcf7_1), 
                                          end(diffloop.mcf7_1), sep="_")

diffloop.mcf7_2 <- makeGRangesFromDataFrame(diffloop.mcf7, seqnames.field="chr_2", start.field="start_2",
                                              end.field="end_2", keep.extra.columns = T)
diffloop.mcf7_2$DiffMC7.anchor_2 <- paste(seqnames(diffloop.mcf7_2), start(diffloop.mcf7_2), 
                                              end(diffloop.mcf7_2), sep="_")
###                                                             ####




##                                    Annotate                                        ##
mcf10a = annotateAnchors2(loops=mcf10a, features=k27.mcf10a, featureName="h3k27peak", featureToAdd = "h3k27peakid", maxgap=0)
mcf10a = annotateAnchors2(loops=mcf10a, features=brca.ccv.gr, featureName="BRCA_CCV", featureToAdd = "ccvid", maxgap=5000)
mcf10a = annotateAnchors2(loops=mcf10a, features=promoters, featureName="HiChIP.Gene", featureToAdd = "gene_name", maxgap=0)
mcf10a = annotateAnchors2(loops=mcf10a, features=brca.twas, featureName="TWAS.Gene", featureToAdd = "Gene_TWAS", maxgap=0)
mcf10a = annotateAnchors2(loops=mcf10a, features=brca.eqtl, featureName="eQTL.Gene", featureToAdd = "eQTL", maxgap=0)
mcf10a = annotateAnchors2(loops=mcf10a, features=diffloop.mcf10a_1, featureName="DiffLoop.MCF10A_1", featureToAdd = "DiffMCF10.anchor_1", maxgap=0)
mcf10a = annotateAnchors2(loops=mcf10a, features=diffloop.mcf10a_2, featureName="DiffLoop.MCF10A_2", featureToAdd = "DiffMCF10.anchor_2", maxgap=0)
mcf10a = annotateAnchorType2(loops=mcf10a, anchor=1) # no ambiguous
mcf10a = annotateAnchorType2(loops=mcf10a, anchor=2) # ambiguous considered E



mcf7 = annotateAnchors2(loops=mcf7, features=k27.mcf7, featureName="h3k27peak", featureToAdd = "h3k27peakid", maxgap=0)
mcf7 = annotateAnchors2(loops=mcf7, features=brca.ccv.gr, featureName="BRCA_CCV", featureToAdd = "ccvid", maxgap=5000)
mcf7 = annotateAnchors2(loops=mcf7, features=promoters, featureName="HiChIP.Gene", featureToAdd = "gene_name", maxgap=0)
mcf7 = annotateAnchors2(loops=mcf7, features=brca.twas, featureName="TWAS.Gene", featureToAdd = "Gene_TWAS", maxgap=0)
mcf7 = annotateAnchors2(loops=mcf7, features=brca.eqtl, featureName="eQTL.Gene", featureToAdd = "eQTL", maxgap=0)
mcf7 = annotateAnchors2(loops=mcf7, features=diffloop.mcf7_1, featureName="DiffLoop.MCF7_1", featureToAdd = "DiffMC7.anchor_1", maxgap=0)
mcf7 = annotateAnchors2(loops=mcf7, features=diffloop.mcf7_2, featureName="DiffLoop.MCF7_2", featureToAdd = "DiffMC7.anchor_2", maxgap=0)

mcf7 = annotateAnchorType2(loops=mcf7, anchor=1) # no ambiguous
mcf7 = annotateAnchorType2(loops=mcf7, anchor=2) # ambiguous considered E
##                                                                                  ##

##                Annotated objects                ##
mcf10a_anno = summary(mcf10a)
mcf10a_anno$loopid = paste(mcf10a_anno$chr_1, mcf10a_anno$start_1, mcf10a_anno$end_1, mcf10a_anno$chr_2, 
                           mcf10a_anno$start_2, mcf10a_anno$end_2, sep="_")

mcf10a_anno$loop.type_details <- apply(mcf10a_anno[,c("anchor.type_1", "anchor.type_2")], 1, 
                                       function(x) paste0(sort(x), collapse = "-"))
mcf10a_anno$loop.type_details_direction = paste(mcf10a_anno$anchor.type_1, mcf10a_anno$anchor.type_2, sep="-")

mcf7_anno = summary(mcf7)
mcf7_anno$loopid = paste(mcf7_anno$chr_1, mcf7_anno$start_1, mcf7_anno$end_1, mcf7_anno$chr_2, 
                           mcf7_anno$start_2, mcf7_anno$end_2, sep="_")

mcf7_anno$loop.type_details <- apply(mcf7_anno[,c("anchor.type_1", "anchor.type_2")], 1, 
                                       function(x) paste0(sort(x), collapse = "-"))
mcf7_anno$loop.type_details_direction = paste(mcf7_anno$anchor.type_1, mcf7_anno$anchor.type_2, sep="-")

mcf10a_anno <- loop_type(mcf10a_anno)
mcf7_anno <- loop_type(mcf7_anno)
## done annotation ##


# 4. is the hichip gene a TWAS gene overlapping anchor 1? # *any* anchor
#MCF10A
t1 = find_match(match1="HiChIP.Gene_1", match2="TWAS.Gene_1", mcf10a_anno, idcol="loopid", split=FALSE, summarizeby = F)
t2 = find_match(match1="HiChIP.Gene_2", match2="TWAS.Gene_2", mcf10a_anno, idcol="loopid", split=FALSE, summarizeby = F)
t3 = find_match(match1="HiChIP.Gene_1", match2="TWAS.Gene_2", mcf10a_anno, idcol="loopid", split=FALSE, summarizeby = F)
t4 = find_match(match1="HiChIP.Gene_2", match2="TWAS.Gene_1", mcf10a_anno, idcol="loopid", split=FALSE, summarizeby = F)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]], t3$loopid[t3[,ncol(t3)]],t4$loopid[t4[,ncol(t4)]]))
mcf10a_anno$TWAS_and_HiChIP = ifelse(mcf10a_anno$loopid %in% id_interest, TRUE, FALSE)

#MCF7
t1 = find_match(match1="HiChIP.Gene_1", match2="TWAS.Gene_1", mcf7_anno, idcol="loopid", split=FALSE, summarizeby = F)
t2 = find_match(match1="HiChIP.Gene_2", match2="TWAS.Gene_2", mcf7_anno, idcol="loopid", split=FALSE, summarizeby = F)
t3 = find_match(match1="HiChIP.Gene_1", match2="TWAS.Gene_2", mcf7_anno, idcol="loopid", split=FALSE, summarizeby = F)
t4 = find_match(match1="HiChIP.Gene_2", match2="TWAS.Gene_1", mcf7_anno, idcol="loopid", split=FALSE, summarizeby = F)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]], t3$loopid[t3[,ncol(t3)]],t4$loopid[t4[,ncol(t4)]]))
mcf7_anno$TWAS_and_HiChIP = ifelse(mcf7_anno$loopid %in% id_interest, TRUE, FALSE)
##                                                                              ##


# is the hichip gene an eQTL gene overlapping anchor 1? # *any* anchor
#MCF10A
t1 = find_match(match1="HiChIP.Gene_1", match2="eQTL.Gene_1", mcf10a_anno, idcol="loopid", split=FALSE, summarizeby = F)
t2 = find_match(match1="HiChIP.Gene_2", match2="eQTL.Gene_2", mcf10a_anno, idcol="loopid", split=FALSE, summarizeby = F)
t3 = find_match(match1="HiChIP.Gene_1", match2="eQTL.Gene_2", mcf10a_anno, idcol="loopid", split=FALSE, summarizeby = F)
t4 = find_match(match1="HiChIP.Gene_2", match2="eQTL.Gene_1", mcf10a_anno, idcol="loopid", split=FALSE, summarizeby = F)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]], t3$loopid[t3[,ncol(t3)]],t4$loopid[t4[,ncol(t4)]]))
mcf10a_anno$eQTL_and_HiChIP = ifelse(mcf10a_anno$loopid %in% id_interest, TRUE, FALSE)

#MCF7
t1 = find_match(match1="HiChIP.Gene_1", match2="eQTL.Gene_1", mcf7_anno, idcol="loopid", split=FALSE, summarizeby = F)
t2 = find_match(match1="HiChIP.Gene_2", match2="eQTL.Gene_2", mcf7_anno, idcol="loopid", split=FALSE, summarizeby = F)
t3 = find_match(match1="HiChIP.Gene_1", match2="eQTL.Gene_2", mcf7_anno, idcol="loopid", split=FALSE, summarizeby = F)
t4 = find_match(match1="HiChIP.Gene_2", match2="eQTL.Gene_1", mcf7_anno, idcol="loopid", split=FALSE, summarizeby = F)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]], t3$loopid[t3[,ncol(t3)]],t4$loopid[t4[,ncol(t4)]]))
mcf7_anno$eQTL_and_HiChIP = ifelse(mcf7_anno$loopid %in% id_interest, TRUE, FALSE)

##       done finding match       ##

mcf7_qtl_twas <- mcf7_anno[mcf7_anno$eQTL_and_HiChIP | mcf7_anno$TWAS_and_HiChIP, ]

## proportion of loops in gene->eqtl or gene->twas 
(dim(mcf7_qtl_twas)[1] / dim(mcf7_anno)[1]) * 100 #0.55%

#table(aux$BRCA_CCV_1 != 'none' & aux$BRCA_CCV_2 != 'none') 
mcf10a_qtl_twas <- mcf10a_anno[mcf10a_anno$eQTL_and_HiChIP | mcf10a_anno$TWAS_and_HiChIP, ]
(dim(mcf10a_qtl_twas)[1] / dim(mcf10a_anno)[1]) * 100 #0.46%

## shared eQTL ##
aux <- OverlapLoop(as.data.frame(mcf7_qtl_twas), as.data.frame(mcf10a_qtl_twas))
shared_eqtl <- mcf7_qtl_twas[unique(c(ab,ba)), ]


### Enhancer only
mcf7_anno_enhancer <- mcf7_anno[mcf7_anno$loop.type == 'E',]
mcf7_anno_enhancer_qtl_twas <- mcf7_anno_enhancer[mcf7_anno_enhancer$eQTL.Gene_1 != 'none' | mcf7_anno_enhancer$eQTL.Gene_2 != 'none', ]
### How many eQTL linked to enhancers???
x <- unique(c(mcf7_anno_enhancer_qtl_twas$eQTL.Gene_1[mcf7_anno_enhancer_qtl_twas$eQTL.Gene_1 != 'none'], 
  mcf7_anno_enhancer_qtl_twas$eQTL.Gene_1[mcf7_anno_enhancer_qtl_twas$eQTL.Gene_2 != 'none'])) 




mcf10a_anno_enhancer <- mcf10a_anno[mcf10a_anno$loop.type == 'E',]
mcf10a_anno_enhancer_qtl_twas <- mcf10a_anno_enhancer[mcf10a_anno_enhancer$eQTL.Gene_1 != 'none' | mcf10a_anno_enhancer$eQTL.Gene_2 != 'none', ]
### How many eQTL linked to enhancers???
y <- unique(c(mcf10a_anno_enhancer_qtl_twas$eQTL.Gene_1[mcf10a_anno_enhancer_qtl_twas$eQTL.Gene_1 != 'none'], 
         mcf10a_anno_enhancer_qtl_twas$eQTL.Gene_1[mcf10a_anno_enhancer_qtl_twas$eQTL.Gene_2 != 'none'])) 

## UNIQUE eQTL GENE In MCF7
x[!x %in% y]
mcf7_anno_enhancer[mcf7_anno_enhancer$eQTL.Gene_1 %in% x[!x %in% y] | mcf7_anno_enhancer$eQTL.Gene_2 %in% x[!x %in% y], ]

## UNIQUE eQTL GENE In MCF10a
y[!y %in% x]



### HICHIP GENE TO eQTL 
mcf7_anno_gene <- mcf7_anno[mcf7_anno$eQTL_and_HiChIP, ]

## EQTL NEGATIVE GENE to CCV
mcf7_anno_qtlneg <- mcf7_anno[mcf7_anno$eQTL_and_HiChIP ==F & (mcf7_anno$BRCA_CCV_1 != 'none' | mcf7_anno$BRCA_CCV_2 != 'none'), ]
mcf7_anno_qtlneg <- mcf7_anno_qtlneg[mcf7_anno_qtlneg$HiChIP.Gene_1 != 'none' | mcf7_anno_qtlneg$HiChIP.Gene_2 != 'none', ]

c(mcf7_anno_qtlneg$HiChIP.Gene_1, mcf7_anno_qtlneg$HiChIP.Gene_2) %>% unique


mcf7_anno_qtlneg_diff <- mcf7_anno_qtlneg[mcf7_anno_qtlneg$DiffLoop.MCF7_1_1 != 'none' | mcf7_anno_qtlneg$DiffLoop.MCF7_1_2 != 'none' |
                   mcf7_anno_qtlneg$DiffLoop.MCF7_2_1 != 'none' | mcf7_anno_qtlneg$DiffLoop.MCF7_2_2 != 'none',]

## DIff looped, Gene->ccv, Enhancer
aux <- mcf7_anno_qtlneg_diff[mcf7_anno_qtlneg_diff$loop.type == 'E-P', ]  

c(aux$HiChIP.Gene_1, aux$HiChIP.Gene_2) %>% unique



mcf10a_anno_qtlneg <- mcf10a_anno[mcf10a_anno$eQTL_and_HiChIP ==F & (mcf10a_anno$BRCA_CCV_1 != 'none' | mcf10a_anno$BRCA_CCV_2 != 'none'), ]
mcf10a_anno_qtlneg <- mcf10a_anno_qtlneg[mcf10a_anno_qtlneg$HiChIP.Gene_1 != 'none' | mcf10a_anno_qtlneg$HiChIP.Gene_2 != 'none', ]

c(mcf10a_anno_qtlneg$HiChIP.Gene_1, mcf10a_anno_qtlneg$HiChIP.Gene_2) %>% unique


mcf10a_anno_qtlneg_diff <- mcf10a_anno_qtlneg[mcf10a_anno_qtlneg$DiffLoop.MCF10A_1_1 != 'none' | mcf10a_anno_qtlneg$DiffLoop.MCF10A_1_2 != 'none' |
                                            mcf10a_anno_qtlneg$DiffLoop.MCF10A_2_1 != 'none' | mcf10a_anno_qtlneg$DiffLoop.MCF10A_2_2 != 'none',]

## DIff looped, Gene->ccv, Enhancer
aux <- mcf10a_anno_qtlneg_diff[mcf10a_anno_qtlneg_diff$loop.type == 'E-P', ]  

c(aux$HiChIP.Gene_1, aux$HiChIP.Gene_2) %>% unique



## save objects
save(mcf7, mcf10a, mcf7_anno, mcf10a_anno, file = '~/Documents/Breast_Hichip/')



##


### Breakdown of number of DiffLoops by loop type
tadiff_by_type.mcf7 <-  table(mcf7_anno$loop.type, (mcf7_anno$DiffLoop.MCF7_1_1 != 'none' | 
                                                      mcf7_anno$DiffLoop.MCF7_1_2 != 'none' | 
                                                      mcf7_anno$DiffLoop.MCF7_2_1 != 'none' | 
                                                      mcf7_anno$DiffLoop.MCF7_2_2 != 'none'))

tadiff_by_type.mcf10a <-  table(mcf10a_anno$loop.type, (mcf10a_anno$DiffLoop.MCF10A_1_1 != 'none' | 
                                                      mcf10a_anno$DiffLoop.MCF10A_1_2 != 'none' | 
                                                      mcf10a_anno$DiffLoop.MCF10A_2_1 != 'none' | 
                                                      mcf10a_anno$DiffLoop.MCF10A_2_2 != 'none'))
###

### Breakdown of number of eQTls by loop type
qtl_by_type.mcf7 <-  table(mcf7_anno$loop.type, (mcf7_anno$eQTL.Gene_1 != 'none' | 
                                                      mcf7_anno$eQTL.Gene_2 != 'none'))

qtl_by_type.mcf10a <-  table(mcf10a_anno$loop.type, (mcf10a_anno$eQTL.Gene_1 != 'none' | 
                                                   mcf10a_anno$eQTL.Gene_2 != 'none'))
###

