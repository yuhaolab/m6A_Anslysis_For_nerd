
# peakGL <-  peakref

.makePeaksummits <- function(peakGL, summit  = 50) {

  peak_name <- peakGL$name
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

  PeakSB <- lapply(PeakL, .makeSBpeak, summit)
  peaklen <- unlist(lapply(PeakSB, length))
  peakcenter <- unlist(PeakSB)
  
  peaksummit <- GRanges(seqnames = rep(peakChr, times = peaklen),
                        IRanges(start = as.numeric(peakcenter), width = 1),
                        strand = rep(peakStrand, times = peaklen))
  peaksummit$name <- rep(peak_name, times = peaklen)

  return(peaksummit)
}

.makeSBpeak <- function(x, extend) {
  x1 <- as.numeric(x[2])
  x2 <- (as.numeric(x[3]) - x1)%/%extend + 1
  x2 <- x1 - 1 + extend*x2
  peakcenter <- (x1):x2
  return(peakcenter)
}


.getpeakmethysummit <- function(peakref, peakrefraw, bams, summit = 50,
                                strandSpecific, ignoreDup, countMultiMappingReads) {
  extend <- summit
  peakbed <- peakref
  bedsaf <- data.frame(GeneID = 1:length(peakbed), Chr = as.character(seqnames(peakbed)), 
                       Start = as.numeric(start(peakbed)), End = as.numeric(end(peakbed)), Strand = as.character(strand(peakbed)))
  x <- featureCounts(bams, 
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
  reads <- x$counts
  reads <- rowSums(reads)
  reads <- pmax(reads, 1)
  reads[countOverlaps(peakref, peakrefraw) == 0] <- 0
  
  peakind <- unique(peakref$name)
  peakpos <- tapply(as.numeric(start(peakref)), peakref$name, c)
  peakreads <- tapply(reads, peakref$name, c)
  peakpos <- peakpos[peakind]
  peakreads <- peakreads[peakind]
  
  centers <- sapply(1:length(peakpos), function(x, peakpos, peakreads, extend, summit) {
    peakpos1 <- peakpos[[x]]
    peakreads1 <- peakreads[[x]]
    len <- length(peakpos1)/extend
    if (len > 2) {
      peakpos1 <- .expendpos(peakpos1, extend)
      peakreads1 <- .expendpos(peakreads1, extend)
      subsmt <- sapply(1:len, function(x, peakpos1, peakreads1) round(weighted.mean(peakpos1[[x]], peakreads1[[x]])), peakpos1, peakreads1)
      
      summitind <- match(subsmt, unlist(peakpos1))
      summitind <- lapply(summitind, function(x, summit) (x-round(summit/2)):(x+round(summit/2)), summit)
      ind <- rep(subsmt, times = unlist(lapply(summitind, length)))
      summitind <- unlist(summitind)
      id <- summitind >=1 & summitind <= len*extend
      summitind <- summitind[id]
      ind <- ind[id]

      summit0 <- tapply(unlist(peakreads1)[summitind], ind, sum)
      summit0 <- names(which.max(summit0))
      
      # summit0 <- lapply(peakreads1, sum)
      # summit0 <- subsmt[which.max(unlist(summit0))]
    } else {
      summit0 <- round(weighted.mean(peakpos1, peakreads1))
    }
    return(as.numeric(summit0))
  }, peakpos, peakreads, extend, summit)
  
  chr <- peakref[match(peakind, peakref$name)]
  chr <- as.character(seqnames(chr))
  
  peaksummit <- GRanges(chr,
                        IRanges(start = centers, width = 1),
                        strand = "*")
  peaksummit$name <- peakind
  # export(peaksummit, "testsummit.bed")
  
  return(peaksummit)
}


# .expendpos <- function(x, extend) {
#   len <- length(x)/extend*2
#   ind <- rep(1:len, each = extend/2)
#   xs <- split(x, ind)
#   xs <- lapply(1:(len - 1), function(i, xs) c(xs[[i]], xs[[i+1]]), xs)
#   return(xs)
# }

.expendpos <- function(x, extend) {
  len <- length(x)/extend
  ind <- rep(1:len, each = extend)
  xs <- split(x, ind)
  return(xs)
}






