

#################################################
#################################################
#################################################
##### function from diffloop
# adaptation: do not collapse anything
loopsMake.mango.null <- function(beddir, samples, mergegap = 500, ext = "all", pad = 0, reduce = TRUE) {
  
  library(foreach)
  library(GenomicRanges)
  library(utils)
  library(dplyr)
  library(readr)
  ext="all"
  loops <- setClass("loops", slots = c(anchors = "GRanges",
                                       interactions = "nim", counts = "nim",
                                       colData = "data.frame", rowData = "data.frame"))
  
  
  ct <- list(col_character(), col_integer(), col_integer(),
             col_character(), col_integer(), col_integer(), col_integer(),
             col_number())
  ct <- list(col_character(), col_integer(), col_integer(),
             col_character(), col_integer(), col_integer(), col_integer(),
             col_number())
  
  # Iterate through files to set up anchors
  anchorsraw <- foreach(sample = samples) %do% {
    fullfile <- file.path(beddir, paste(sample, "interactions", ext, "mango", sep = "."))
    bt <- read_delim(fullfile, "\t", col_types = ct, col_names = FALSE)
    plyr::rbind.fill(bt[, 1:3], setNames(bt[, 4:6], names(bt[, 1:3])))
  }
  
  if (reduce) {
    anchors <- reduce(padGRanges(makeGRangesFromDataFrame(do.call(rbind,
                                                                  anchorsraw), ignore.strand = TRUE, seqnames.field = "X1",
                                                          start.field = "X2", end.field = "X3"), pad = pad), min.gapwidth = mergegap)
  }
  if (!reduce) {
    anchors <- padGRanges(makeGRangesFromDataFrame(do.call(rbind,
                                                           anchorsraw), ignore.strand = TRUE, seqnames.field = "X1",
                                                   start.field = "X2", end.field = "X3"), pad = pad)
  }
  
  
  petlist <- foreach(sample = samples) %do% {
    fullfile <- file.path(beddir, paste(sample, "interactions", ext, "mango", sep = "."))
    bt <- read_delim(fullfile, "\t", col_types = ct, col_names = FALSE)
    
    anchors_left = makeGRangesFromDataFrame(bt[, 1:3], ignore.strand = TRUE, seqnames.field = "X1",
                                            start.field = "X2", end.field = "X3")
    anchors_right = makeGRangesFromDataFrame(bt[, 4:6], ignore.strand = TRUE, seqnames.field = "X4",
                                             start.field = "X5", end.field = "X6")
    counts = bt[, 7]
    
    leftanchor = 1:(length(anchors)/2)
    rightanchor = ((length(anchors)/2)+1):length(anchors)
    df <- data.frame(left = leftanchor, right = rightanchor)
    g <- as.data.frame(dplyr::group_by(df, left, right))
    d <- cbind(g, counts)
    colnames(d) <- c("left", "right", "counts")
    dag <- d
    #dag <- aggregate(counts ~ left + right, FUN = sum, data=d)
    colnames(dag) <- c("left", "right", sample)
    dag
  }
  
  .full_join <- function(a, b) {
    as.data.frame(dplyr::full_join(a, b, by = c("left", "right")))
  }
  
  # Map Counts
  pets <- Reduce(.full_join, petlist)
  iraw <- pets[, c("left", "right")]
  iraw <- t(apply(iraw, 1, function(x) {
    if (x[1] < x[2]) {
      x
    } else {
      x[c(2, 1)]
    }
  }))
  interactions <- iraw[order(iraw[, 1], iraw[, 2]), ]
  colnames(interactions) <- c("left", "right")
  counts <- as.matrix(pets[, -c(1:2)])[order(iraw[, 1], iraw[, 2]), ]
  counts[is.na(counts)] <- 0
  counts <- as.matrix(counts, ncol = length(samples))
  colnames(counts) <- samples
  
  # Initialize rowData slot (with loop widths)
  w <- (start(anchors[interactions[, 2]]) + end(anchors[interactions[, 2]]))/2 -
    (start(anchors[interactions[, 1]]) + end(anchors[interactions[, 1]]))/2
  w[w < 0] <- 0
  rowData <- as.data.frame(as.integer(w))
  colnames(rowData) <- c("loopWidth")
  
  # Remove 'chr' from anchors
  seqlevels(anchors) <- gsub("^chr(.*)$", "\\1", seqlevels(anchors))
  
  # Remove rownames from matrices
  row.names(interactions) <- NULL
  row.names(counts) <- NULL
  
  # Initialize colData slot
  groups <- rep("group1", length(samples))
  if(length(samples) == 1){
    sizeFactor <- 1
  } else {
    lc <- log2(counts)
    keep <- rowSums(counts > 0) == ncol(lc)
    lc <- lc[keep, ]
    target <- 2^rowMeans(lc)
    sizeFactor <- colMedians(sweep(2^lc, 1, target, FUN = "/"), na.rm = TRUE)
  }
  dfcd <- data.frame(sizeFactor, groups)
  rownames(dfcd) <- samples
  
  # Create loops object
  dlo <- loops()
  slot(dlo, "anchors", check = TRUE) <- anchors
  slot(dlo, "interactions", check = TRUE) <- interactions
  slot(dlo, "counts", check = TRUE) <- counts
  slot(dlo, "colData", check = TRUE) <- dfcd
  slot(dlo, "rowData", check = TRUE) <- rowData
  
  return(dlo)
}


annotateAnchors2 <- function(loops=all, features=eqtl_sig_gr, featureName="EQTL.TCGA", featureToAdd = "eQTLs", maxgap=0) {
  if (!(featureToAdd %in% names(mcols(features)))) stop("Do not have featureToAdd in data")
  lto = loops
  lto.df <- summary(lto)
  Lanchors <- GRanges(lto.df$chr_1, IRanges(lto.df$start_1, lto.df$end_1))
  Ranchors <- GRanges(lto.df$chr_2, IRanges(lto.df$start_2, lto.df$end_2))
  
  # Determine if right anchor overlaps GWAS region
  Rhits.p <- suppressWarnings(findOverlaps(features, Ranchors,
                                           maxgap = maxgap))
  Rvalues.p <- rep(FALSE, dim(lto.df)[1])
  Rvalues.p[unique(subjectHits(Rhits.p))] <- TRUE
  
  # Determine if left anchor overlaps GWAS region
  Lhits.p <- suppressWarnings(findOverlaps(features, Lanchors,
                                           maxgap = maxgap))
  Lvalues.p <- rep(FALSE, dim(lto.df)[1])
  Lvalues.p[unique(subjectHits(Lhits.p))] <- TRUE
  
  # Aggregate TSS
  Rtss <- data.frame(Rhits.p)
  Rtss <- cbind(Rtss, mcols(features[Rtss$queryHits]))
  Ltss <- data.frame(Lhits.p)
  Ltss <- cbind(Ltss, mcols(features[Ltss$queryHits]))
  # collapse
  r_ttss <- unique(Rtss[,c(2,which(names(Rtss) == featureToAdd))])
  r_tss <- aggregate(r_ttss[,featureToAdd]~subjectHits,paste,collapse=",",data=r_ttss)
  names(r_tss)[2] = featureToAdd
  
  l_ttss <- unique(Ltss[,c(2,which(names(Ltss) == featureToAdd))])
  l_tss <- aggregate(l_ttss[,featureToAdd]~subjectHits,paste,collapse=",",data=l_ttss)
  names(l_tss)[2] = featureToAdd
  
  anchor.1 <- rep("none", dim(lto)[2]) # record left is 1
  anchor.2 <- rep("none", dim(lto)[2]) # record right as 2
  
  anchor.1[l_tss$subjectHits] <- l_tss[,featureToAdd]
  anchor.2[r_tss$subjectHits] <- r_tss[,featureToAdd]
  
  lto@rowData$anchor.1 <- anchor.1
  lto@rowData$anchor.2 <- anchor.2
  
  colnames(lto@rowData)[which(colnames(lto@rowData) == "anchor.1")] = paste(featureName, "_1", sep="")
  colnames(lto@rowData)[which(colnames(lto@rowData) == "anchor.2")] = paste(featureName, "_2", sep="")
  
  return(lto)
}


annotateAnchorType2 <- function(
  loops,
  h3k27peakid = "h3k27peak",
  gene = "HiChIP.Gene",
  anchor=1) {
  
  df_anno <- summary(loops)
  loop.types <- 1:nrow(df_anno)
  description <- rep("none", length(loop.types))
  
  h3k27peakid = paste(h3k27peakid, anchor, sep="_")
  gene = paste(gene, anchor, sep="_")
  
  if ( !h3k27peakid %in% names(df_anno) ) stop("h3k27peakid column not in data")
  if ( !gene %in% names(df_anno) ) stop("gene column not in data")
  
  p1 = (
    df_anno[,h3k27peakid] != "none" &
      df_anno[,gene] !="none"
  )
  p2 = (
    df_anno[,h3k27peakid] == "none" &
      df_anno[,gene] !="none"
  )
  if ( length(intersect(which(p1),which(p2))) > 0) stop("Promoter definitions overlapping")
  p = p1 | p2
  #length(which(p))
  
  e = (
    df_anno[,h3k27peakid] != "none" &
      df_anno[,gene] =="none"
  )
  
  #if ( length(intersect(which(e1),which(e2)))  > 0) stop("Enhancer definitions overlapping")
  #e = e1 | e2
  #length(which(e))
  
  if ( length(intersect(which(p),which(e))) > 0) stop()
  
  other = (
    df_anno[,h3k27peakid] == "none" &
      df_anno[,gene] =="none"
  )
  
  if ( length(intersect(which(p),which(other))) > 0) stop("other definitions overlapping")
  if ( length(intersect(which(e),which(other))) > 0) stop("other definitions overlapping")
  
  
  # sum all: should be equal to number of loops
  description[p] <- "P"
  description[e] <- "E"
  description[other] <- "O"
  
  if ( "none" %in% names(table(description)) )
    print("Unclassified cases (none): look what these are!")
  loops@rowData$type <- description
  names(loops@rowData)[ncol(loops@rowData)] = paste("anchor.type_", anchor, sep="")
  #names(loops@rowData)[ncol(loops@rowData)] = paste("anchor.type.noamb_", anchor, sep="")
  return(loops)
  
}

# if summarizeby, return other columns grouped by loopid or anchorid
# if not summarizeby, return the original dataset with added true/false for matching
find_match = function(match1="HiChIP.Gene_1", match2="Gene.TWAS.anchor_1", df, idcol, split=FALSE, split_match = 2, summarizeby = TRUE, equal = TRUE) {
  if (!(match1 %in% colnames(df) )) stop("Column ", match1, " is missing")
  if (!(match2 %in% colnames(df) )) stop("Column ", match2, " is missing")
  if (!(idcol %in% colnames(df) )) stop("Column ", idcol, " is missing")
  if (any(grepl(",", df[,idcol]))) stop("idcol has already been grouped")
  
  
  t = df
  t = t[t[,match1] !="none" | t[,match2]!="none",] # to help with speed
  t = separate_rows(t, match1, sep = ",")
  t = separate_rows(t, match2, sep = ",")
  if (split & split_match == 1) {
    t[,match1] = as.character(lapply(strsplit(as.character(t[,match1]), ":", fixed=TRUE), "[", 2))
  }
  if (split & split_match == 2) {
    t[,match2] = as.character(lapply(strsplit(as.character(t[,match2]), ":", fixed=TRUE), "[", 2))
  }
  
  t = t[t[,match1] !="none",]
  t = t[t[,match2] !="none",]
  t = t[!is.na(t[,match1]),]
  t = t[!is.na(t[,match1]),]
  
  if (equal) {
    t =t[t[,match1] == t[,match2],]
  } else {
    t =t[t[,match1] != t[,match2],]
  }
  t = unique(t)
  if (!summarizeby) {
    id_interest = as.data.frame(unique(t[,idcol]))
    df$temp = ifelse(df[,idcol] %in% id_interest[, 1], TRUE, FALSE)
    names(df)[ncol(df)] = paste(match1, match2, sep="_")
    return(df)
  } else {
    # collapse by id
    # print(cbind(t[,match1], t[,match2]))
    library(plyr)
    # tss <- aggregate(.~t[,idcol], data = t, function(x) paste0(unique(x), collapse=","))
    # https://stackoverflow.com/questions/53341434/group-by-a-column-and-collapse-all-other-columns-without-na
    tss = data.frame (t %>%
                        mutate_if(is.factor, as.character)  %>%
                        mutate_all(~replace(., .=='NA', NA)) %>%
                        group_by(t[,idcol]) %>%
                        summarize_all(~paste(unique(na.omit(.)), collapse = ',')) )
    #summarise_all(funs(unique(paste(ifelse(is.na(.), "null", .), collapse = ",")))) )
    return(tss)
  }
}

loop_type <- function(df_anno, with_ambiguous = F) {
  if (with_ambiguous) {
    # if with_ambiguous = TRUE, loop.type_details is 9 categories, now collapse into 4 categories
    e.loops = df_anno$loop.type_details %in% c("A-E", "E-E", "E-O")
    p.loops = df_anno$loop.type_details %in% c("A-P", "O-P", "P-P")
    ep.loops = df_anno$loop.type_details %in% c("E-P")
    a.loops = df_anno$loop.type_details %in% c("A-A", "A-O")
    df_anno$loop.type = "none"
    df_anno$loop.type[e.loops] <- "E"
    df_anno$loop.type[p.loops] <- "P"
    df_anno$loop.type[ep.loops] <- "E-P"
    df_anno$loop.type[a.loops] <- "A"
  }
  if (!with_ambiguous) {
    e.loops = df_anno$loop.type_details %in% c("E-E", "E-O")
    p.loops = df_anno$loop.type_details %in% c("O-P", "P-P")
    ep.loops = df_anno$loop.type_details %in% c("E-P")
    df_anno$loop.type = "none"
    df_anno$loop.type[e.loops] <- "E"
    df_anno$loop.type[p.loops] <- "P"
    df_anno$loop.type[ep.loops] <- "E-P"
  }
  df_anno$loop.type_details = gsub("O-P", "P-O", df_anno$loop.type_details)
  return(df_anno)
  
}