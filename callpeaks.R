source("/media/yh/Nano2/YL/m6a/peakcallingmodule/modulesFun.R")
library(GenomicFeatures)
library(rtracklayer)
library(Rsubread)
library(exomePeak)

txdb <- makeTxDbFromGFF("/media/yh/Nano2/refnopt/Arabidopsis_thaliana.TAIR10.50.gtf")
genegr <- genes(txdb)

# setwd("/media/yh/Nano2/YL/m6a/peakcalling/peak50bs")
# peakpoint <- .makePeaksummits(genegr, extend = 50, steps = 50)

# setwd("/media/yh/Nano2/YL/m6a/peakcalling/peak25bs")
# peakpoint <- .makePeaksummits(genegr, extend = 25, steps = 25)

# setwd("/media/yh/Nano2/YL/m6a/peakcalling/peak10b1s")
# library(data.table)
# peakpoint <- .makePeaksummits(genegr, extend = 10, steps = 1)
# peakpoint <- resize(peakpoint, width = 1, fix = "start", ignore.strand = T)
# peakreads <- fread("peak_counts.txt")
# peakreads <- data.frame(peakreads)

# setwd("/media/yh/Nano2/YL/m6a/peakcalling/peak100b50s")
# peakpoint <- .makePeaksummits(genegr, extend = 100, steps = 50)
# peakpoint <- resize(peakpoint, width = 50, fix = "start", ignore.strand = T)

# setwd("/media/yh/Nano2/YL/m6a/peakcalling/peak100b25s/")
# peakpoint <- .makePeaksummits(genegr, extend = 100, steps = 25)
# peakpoint <- resize(peakpoint, width = 25, fix = "start", ignore.strand = T)

# setwd("/media/yh/Nano2/YL/m6a/peakcalling/peak50b20s")
# peakpoint <- .makePeaksummits(genegr, extend = 50, steps = 20)
# peakpoint <- resize(peakpoint, width = 20, fix = "start", ignore.strand = T)

setwd("/media/yh/Nano2/YL/m6a/peakcalling/peak50b25s")
peakpoint <- .makePeaksummits(genegr, extend = 50, steps = 25)
peakpoint <- resize(peakpoint, width = 25, fix = "start", ignore.strand = T)

peakreads <- read.delim("peak_counts.txt")
genereads <- read.delim("gene_counts.txt")

## get total reads
# sf <- colSums(genereads[,-(1:6)])
# wtiptotal <- sum(sf[grep("WT_MeRIP",names(sf))])
# wtinputtotal <- sum(sf[grep("WT_input",names(sf))])
# mutiptotal <- sum(sf[grep("Mu_MeRIP",names(sf))])
# mutinputtotal <- sum(sf[grep("Mu_input",names(sf))])

sf <- read.delim("gene_counts_summary.txt")
sf <- sf[1,-1]
wtiptotal <- sum(sf[grep("WT_MeRIP",names(sf))])
wtinputtotal <- sum(sf[grep("WT_input",names(sf))])
mutiptotal <- sum(sf[grep("Mu_MeRIP",names(sf))])
mutinputtotal <- sum(sf[grep("Mu_input",names(sf))])


fcpre <- 2
# fcpre <- 1+1/log(candipeak$ip, base = exp(1))

## peak calling for WT
candipeak <- peakreads[,1:6]
candipeak$ind <- 1:nrow(candipeak)
candipeak$ip <- rowSums(peakreads[,grep("WT_MeRIP", names(peakreads))])
candipeak$input <- rowSums(peakreads[,grep("WT_input", names(peakreads))])

## test
re <- ctest(candipeak$ip, candipeak$input, wtiptotal, wtinputtotal, minimal_counts_in_fdr = 10)
re <- data.frame(fc=exp(re$log.fc), p=exp(re$log.p), fdr=exp(re$log.fdr))
candipeak <- cbind(candipeak, re)

## significant
sigid <- candipeak$fdr < 0.05 & candipeak$fc >= fcpre
peakwt <- candipeak[sigid,]

length(unique(peakwt$GeneID))


## mutant
candipeak <- peakreads[,1:6]
candipeak$ind <- 1:nrow(candipeak)
candipeak$ip <- rowSums(peakreads[,grep("Mu_MeRIP", names(peakreads))])
candipeak$input <- rowSums(peakreads[,grep("Mu_input", names(peakreads))])

## test
re <- ctest(candipeak$ip, candipeak$input, mutiptotal, mutinputtotal, minimal_counts_in_fdr = 10)
re <- data.frame(fc=exp(re$log.fc), p=exp(re$log.p), fdr=exp(re$log.fdr))
candipeak <- cbind(candipeak, re)

## sig
sigid <- candipeak$fdr < 0.05 & candipeak$fc >= fcpre
peakmut <- candipeak[sigid,]

length(unique(peakmut$GeneID))


## test flc
flcwt <- peakwt[peakwt$GeneID == "AT5G10140",]
flcmut <- peakmut[peakmut$GeneID == "AT5G10140",]


peakmac <- import("/media/yh/Nano2/YL/m6a/old/peakcalling/macspeak_comb/WT_peaks.narrowPeak")
peak0 <- makeGRangesFromDataFrame(peakwt)
sum(countOverlaps(peakmac, peak0) != 0)/length(peakmac)


## -----------------------------------------------------------------------------------------------------------------------------------------
## output
grwt <- .writepeaks(peakwt, peakpoint, 100, 2, "wt")
grmut <- .writepeaks(peakmut, peakpoint, 100, 2, "mut")

summary(grwt$ip)
summary(grmut$ip)
length(unique(grwt$GeneID))
length(unique(grmut$GeneID))

sum(countOverlaps(grwt, grmut) != 0)

library(ggvenn)

a <- list(`wt` = unique(grwt$GeneID),
          `mut` = unique(grmut$GeneID))
ggvenn(a, auto_scale = T)

length(grwt)
length(grmut)
sum(countOverlaps(peakmac, grwt) != 0)/length(peakmac)


peak <- mcols(grwt)
peak <- data.frame(peak)













