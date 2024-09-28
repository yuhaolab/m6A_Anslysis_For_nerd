library(DESeq2)
library(ggplot2)


# gtf = "/home/yh/3-tool-ref/ref/Arabidopsis_thaliana.TAIR10.50.gtf"
# recenter = T
# summit = 30
# extend = 100
# countgene = T
# countpeak = T
# normalize = "gene_total" ## other: "gene_deseq2", "peak_total", "peak_deseq2"
# diffmethod = "exomepeak2" ## other: exomepeak
# outdir = NA
# strandSpecific = 0
# ignoreDup = F
# countMultiMappingReads = F
# savename = "tr_wt"
# fc_thr = 0.5
# q_thr = 0.05
# p_thr = NA
# centerind = "methy" ## other: "p", "fc"


diffpeaksummit <- function(peakwt, peakmut, wt_ip, wt_input, tr_ip, tr_input,
                           gtf = "/home/yh/3-tool-ref/ref/Arabidopsis_thaliana.TAIR10.50.gtf",
                           recenter = T,
                           summit = 30, 
                           extend = 100,
                           countgene = T,
                           countpeak = T,
                           normalize = "gene_total",
                           diffmethod = "exomepeak2",
                           outdir = NA,
                           strandSpecific = 0,
                           ignoreDup = F,
                           countMultiMappingReads = F,
                           savename = "tr_wt",
                           fc_thr = 1,
                           q_thr = 0.05,
                           p_thr = NA,
                           centerind = "methy") {
  
  if (is.na(outdir)) {outdir <- getwd()}
  
  peakref <- reduce(c(peakmut, peakwt))
  peakref$name <- paste0("p_", 1:length(peakref))
  
  if (recenter) {
    peakrefraw <- peakref
    peakref <- .makePeaksummits(peakref, summit)
    
    if (centerind == "methy") {
      bams <- c(wt_ip, tr_ip)
      peakref <- .getpeakmethysummit(peakref, peakrefraw, bams, summit, strandSpecific, ignoreDup, countMultiMappingReads)
    }
    peaksummit <- as.numeric(start(peakref))
    peakref <- resize(peakref, width = extend + 1, fix = "center")
    # export(peakref, "testextend.bed")
  }

  if (countpeak) {
    .getpeakreads(peakref = peakref, wt_ip = wt_ip, wt_input = wt_input, tr_ip = tr_ip, tr_input = tr_input,
                  gtf = gtf, countgene = countgene, outdir = outdir,
                  strandSpecific = strandSpecific, ignoreDup = ignoreDup, countMultiMappingReads = countMultiMappingReads)
  }
  
  
  ## get diff --------------------------------------------------------------------------------------------------------------------------------------
  diffpeak <- peakref
  reads <- read.table(paste0(outdir, "/peak_counts.txt"), header = T, stringsAsFactors = F)
  peakname <- reads$GeneID
  rownames(reads) <- make.unique(reads$GeneID)
  reads <- reads[,7:ncol(reads)]
  
  cond = rep(c("wt_ip", "wt_input", "tr_ip", "tr_input"), times = c(length(wt_ip), length(wt_input), length(tr_ip), length(tr_input)))
  
  if (grepl("gene", normalize)) {
    genereads <- read.table(paste0(outdir, "/gene_counts.txt"), header = T, stringsAsFactors = F)
    rownames(genereads) = genereads$GeneID
    genereads = genereads[,7:ncol(genereads)]
    
    if (grepl("total", normalize)) {
      print("normalized with total number of gene reads......")
      sf <- colSums(genereads)
      # sf <- read.table(paste0(outdir, "/gene_counts_summary.txt"), header = T, stringsAsFactors = F)
      # sf <- as.numeric(sf[1,-1])
      sf <- sf/median(sf)
    } else {
      print("normalized with gene reads using DESeq2......")
      sf <- estimateSizeFactorsForMatrix(as.matrix(genereads))
    }
  } else if (normalize == "peak_total") {
    print("normalized with total number of peak reads......")
    sf <- colSums(reads)
    sf <- sf/median(sf)
  } else {
    print("normalized with peak reads using DESeq2......")
    sf <- estimateSizeFactorsForMatrix(as.matrix(reads))
  }
  
  if (diffmethod == "exomepeak2") { re <- .getdiff2(reads, cond, diffpeak, sf) }
  
  hist(re$log2fc)
  
  if (!is.na(p_thr)) {
    print("filter significant diff peak with pvalue and log2 foldchange......")
    treid <- re$p < p_thr & re$log2fc > fc_thr
    uneid <- re$p < p_thr & re$log2fc < -fc_thr
  } else {
    print("filter significant diff peak with qvalue and log2 foldchange......")
    treid <- re$fdr < q_thr & re$log2fc > fc_thr
    uneid <- re$fdr < q_thr & re$log2fc < -fc_thr
  }
  
  re$WT_enhance <- uneid
  re$mut_enhance <- treid
  
  wtunique <- countOverlaps(diffpeak, peakwt) != 0 & countOverlaps(diffpeak, peakmut) == 0
  mutunique <- countOverlaps(diffpeak, peakmut) != 0 & countOverlaps(diffpeak, peakwt) == 0
  common <- countOverlaps(diffpeak, peakmut) != 0 & countOverlaps(diffpeak, peakwt) != 0
  
  re$WT_unique <- wtunique
  re$mut_unique <- mutunique
  re$common <- common
  
  ## get methy level
  rfc <- t(t(reads)/sf)
  colnames(rfc) <- paste(cond, c(1:length(wt_ip), 1:length(wt_input), 1:length(tr_ip), 1:length(tr_input)), sep = ".")
  wtmethy <- rowMeans(rfc[,grep("wt_ip", colnames(rfc))])/rowMeans(rfc[,grep("wt_input", colnames(rfc))])
  mutmethy <- rowMeans(rfc[,grep("tr_ip", colnames(rfc))])/rowMeans(rfc[,grep("tr_input", colnames(rfc))])
  re <- cbind(re, rfc, wtmethy, mutmethy)
  
  
  if (recenter) {
    
    # if (centerind == "p") {
    #   print("re-centered peak based on diff pvalue......")
    #   indp <- re$p
    #   ind <- tapply(indp, peakname, min)
    # } else {
    #   if (centerind == "fc") {
    #     print("re-centered peak based on fold change......")
    #     indp <- abs(re$log2fc)
    #   }
    #   if (centerind == "methy") {
    #     print("re-centered peak based on methy level......")
    #     indp <- rowMeans(rfc[,grep("_ip", colnames(rfc))])
    #   }
    #   ind <- tapply(indp, peakname, max)
    # }
    # 
    # ind <- ind[unique(peakname)]
    # ind <- match(paste(names(ind), as.numeric(ind)), paste(peakname, indp))
    # 
    # export(resize(peakref[ind], width = 1, fix = "center"),  paste0(savename, "_summit.bed"))
    # 
    # re$summit <- peaksummit
    # re <-re[ind,]
    
    peaksummit <- resize(peakref, width = 1, fix = "center")
    export(peaksummit,  paste0(savename, "_summit.bed"))
    
    re$summit <- as.numeric(start(peaksummit))
    diffpeak <- peakrefraw
  }
  
  xls <- data.frame(
    chr = as.character(seqnames(diffpeak)),
    start = start(diffpeak),
    end = end(diffpeak),
    strand = as.character(strand(diffpeak)),
    name = diffpeak$name
  )
  
  ## anno gene
  txdb <- makeTxDbFromGFF(gtf)
  
  genecomp <- genes(txdb)
  ind <- findOverlaps(diffpeak, genecomp)
  peakgene <- tapply(genecomp$gene_id[subjectHits(ind)], queryHits(ind), unique)
  peakgene <- lapply(peakgene, paste, collapse = ",")
  
  xls$gene <- NA
  xls$gene[as.numeric(names(peakgene))] <- unlist(peakgene)
  xls <- cbind(xls, re)
  
  sum(xls$mut_enhance)
  sum(xls$WT_enhance)
  sum(xls$p < 0.05)/nrow(xls)
  length(unique(xls[xls$mut_enhance | xls$WT_enhance,]$name))
  
  # sigid <- uneid | treid
  write.table(xls, file = paste0(savename, "_diff.xls"), sep = "\t", quote = F, col.names = T, row.names = F)
  export(diffpeak,  paste0(savename, "_diff.bed"))
  # write.table(xls[sigid,], file = paste0(savename, "_sigdiff.xls"), sep = "\t", quote = F, col.names = T, row.names = F)
  
  return(xls)
}


## get diff exomepeak2 ##################################################################################################################
.getdiff2 <- function(reads, cond, diffpeak, sf = NULL) {
  
  ## make condition
  cond1 <- sub(".*_", "", cond)
  cond2 <- sub("_.*", "", cond)
  se <- data.frame(IP_input = cond1, Perturbation = cond2)
  se$IP_input <- relevel(factor(se$IP_input),"input")
  se$Perturbation <- relevel(factor(se$Perturbation),"wt")
  
  ## make deseq2 dataset
  dds <- DESeqDataSetFromMatrix(reads, se, ~ IP_input * Perturbation)
  if (is.null(sf)) {
    dds <- estimateSizeFactors(dds)
  } else {
    dds$sizeFactor <- sf
  }
  
  ## poison distribution
  dispersions(dds) <- 0
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  # plotMA(res, alpha = 0.05)
  
  pvalue <- res$pvalue
  pvalue[is.na(pvalue)] <- 1
  log2fc <- res$log2FoldChange
  re <- data.frame(reads=res$baseMean, log2fc = log2fc, p = pvalue, fdr = p.adjust(pvalue, method = "BH"))
  
  return(re)
}












