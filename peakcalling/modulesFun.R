# peakGL <- genegr[1:10]
# extend  = 50
# steps = 50
# peakbed <- peakpoint


.makePeaksummits <- function(peakGL, extend  = 50, steps = 50) {
  
  peak_name <- peakGL$gene_id
  peak_name[is.na(peak_name)] <- "Unk"
  
  peakStart <- as.numeric(start(peakGL))
  peakEnd <- as.numeric(end(peakGL))
  peakChr <- as.character(seqnames(peakGL))
  peakStrand <- as.character(strand(peakGL))
  
  len = length(peakGL)
  PeakM <- rbind(peakChr, peakStart, peakEnd, peakStrand)
  PeakL <- as.vector(PeakM)
  ind <- rep(1:len, each = 4)
  PeakL <- split(PeakL, ind)
  names(PeakL) <- peak_name
  
  PeakSB <- lapply(PeakL, .makeSBpeak, steps)
  peaklen <- unlist(lapply(PeakSB, length))
  peakcenter <- unlist(PeakSB)
  
  peaksummit <- GRanges(seqnames = rep(peakChr, times = peaklen),
                        IRanges(start = as.numeric(peakcenter), width = extend),
                        strand = rep(peakStrand, times = peaklen))
  peaksummit$name <- rep(peak_name, times = peaklen)
  
  return(peaksummit)
}

.makeSBpeak <- function(x, steps) {
  x1 <- as.numeric(x[2])
  x2 <- as.numeric(x[3])
  peakcenter <- seq(x1, x2, by = steps)
  return(peakcenter)
}


.getpeakreads <- function(peakbed, ipbam, inputbam,
                          gtf = "/media/yh/Nano2/refnopt/Arabidopsis_thaliana.TAIR10.50.gtf",
                          countgene = T,
                          outdir = NA,
                          strandSpecific = 0,
                          ignoreDup = F,
                          countMultiMappingReads = F) {
  
  if (is.na(outdir)) {outdir <- getwd()}
  
  ## count ---------------------------------------------------------------------------------------------------------------------------------------
  
  files <- c(ipbam, inputbam)
  bedsaf <- data.frame(GeneID = 1:length(peakbed), Chr = as.character(seqnames(peakbed)), 
                       Start = as.numeric(start(peakbed)), End = as.numeric(end(peakbed)), Strand = as.character(strand(peakbed)))
  x <- featureCounts(files, 
                     annot.ext = bedsaf,
                     isGTFAnnotationFile = F,
                     countMultiMappingReads = countMultiMappingReads,
                     allowMultiOverlap = T,
                     isPairedEnd = T,
                     requireBothEndsMapped = T,
                     countChimericFragments = F,
                     ignoreDup = ignoreDup,
                     strandSpecific = strandSpecific,
                     nthreads = 32)
  
  res <- cbind(x$annotation, x$counts)
  names(res)[7:ncol(res)] <- sub(".markdup.*", "", names(res)[7:ncol(res)])
  ress <- x$stat
  names(ress)[2:ncol(ress)] <- sub(".markdup.*", "", names(res)[7:ncol(res)])
  res$GeneID <- peakbed$name
  
  write.table(res, file = paste0(outdir, "/peak_counts.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(ress, file = paste0(outdir, "/peak_counts_summary.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  
  if (countgene) {
    ## count gene
    x <- featureCounts(files, 
                       annot.ext = gtf,
                       isGTFAnnotationFile = T,
                       GTF.featureType = "gene",
                       countMultiMappingReads = countMultiMappingReads,
                       allowMultiOverlap = T,
                       isPairedEnd = T,
                       requireBothEndsMapped = T,
                       countChimericFragments = F,
                       ignoreDup = countMultiMappingReads,
                       strandSpecific = strandSpecific,
                       nthreads = 32)
    
    res <- cbind(x$annotation, x$counts)
    names(res)[7:ncol(res)] <- sub(".markdup.*", "", names(res)[7:ncol(res)])
    ress <- x$stat
    names(ress)[2:ncol(ress)] <- sub(".markdup.*", "", names(res)[7:ncol(res)])
    
    write.table(res, file = paste0(outdir, "/gene_counts.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(ress, file = paste0(outdir, "/gene_counts_summary.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  }
}


## outpot ---------------------------------------------------------------------------------------------------------------------------------
.writepeaks <- function(peak, peakpoint, minpeakwidth = 100, minfc = 3, outname = "peak") {
  gr <- peakpoint[peak$ind]
  gr <- reduce(gr)
  gr <- gr[width(gr) >= minpeakwidth]
  
  # sum(countOverlaps(peakmac, gr) != 0)/length(peakmac)
  
  grsig <- makeGRangesFromDataFrame(peak)
  ind <- findOverlaps(gr, grsig)
  
  fc <- tapply(peak$fc[subjectHits(ind)], queryHits(ind), max)
  id <- match(paste(names(fc), fc), paste(queryHits(ind), peak$fc[subjectHits(ind)]))
  id <- subjectHits(ind)[id]
  peaksig <- peak[id,]
  
  names(peaksig)[2:5] <- paste0("thin_", names(peaksig)[2:5])
  ind <- tapply(peak$ind[subjectHits(ind)], queryHits(ind), paste, collapse = ",")
  peaksig$Allind <- unlist(ind)
  mcols(gr) <- peaksig
  
  gr <- gr[gr$fc >= minfc]
  # sum(countOverlaps(peakmac, gr) != 0)/length(peakmac)
  
  # export(gr, "wtpeak.bed")
  # grwt <- gr
  export(gr, paste0(outname, "_peak.bed"))
  # grmut <- gr
  
  xls <- mcols(gr)
  xls <- data.frame(xls)
  names(xls)[2:5] <- sub("thin_", "", names(xls)[2:5])
  xls$Start <- start(gr)
  xls$End <- end(gr)
  xls$Length <- width(gr)
  
  write.table(xls, file = paste0(outname, "_peak.xls"), quote = F, sep = "\t", col.names = T, row.names = F)
  
  return(gr)
}








