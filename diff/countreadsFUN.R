library(GenomicFeatures)
library(rtracklayer)
library(Rsubread)
# cond = c("wt_ip", "wt_ip", "wt_input", "wt_input", "tr_ip", "tr_ip", "tr_input", "tr_input")

.getpeakreads <- function(peakref, wt_ip, wt_input, tr_ip, tr_input,
                          gtf = "/home/yh/3-tool-ref/ref/Arabidopsis_thaliana.TAIR10.50.gtf",
                          countgene = T,
                          outdir = NA,
                          strandSpecific = 0,
                          ignoreDup = F,
                          countMultiMappingReads = F) {
  
  if (is.na(outdir)) {outdir <- getwd()}
  
  ## count ---------------------------------------------------------------------------------------------------------------------------------------
  peakbed <- peakref
  
  files <- c(wt_ip, wt_input, tr_ip, tr_input)
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