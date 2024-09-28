source("/media/yh/Nano2/YL/m6a/diffpeakmodule/countreadsFUN.R")
source("/media/yh/Nano2/YL/m6a/diffpeakmodule/cutpeaksFUN.R")
source("/media/yh/Nano2/YL/m6a/diffpeakmodule/diffpeaksummitFUN.R")

setwd("/media/yh/Nano2/YL/m6a/bamnpt")
bamfile <- list.files(path = getwd(), pattern = "[.]bam$", full.names = T)
inputbam <- bamfile[grep("input", bamfile)]
ipbam <- bamfile[grep("MeRIP", bamfile)]

wt_ip <- ipbam[grep("WT_", ipbam)]
wt_input <- inputbam[grep("WT_", inputbam)]
tr_ip <- ipbam[grep("Mu_", ipbam)]
tr_input <- inputbam[grep("Mu_", inputbam)]

setwd("/media/yh/Nano2/YL/m6a/diffpeak")
peakmut <- import("../peakcalling/peak50b25s/mut_peak.bed")
peakwt <- import("../peakcalling/peak50b25s/wt_peak.bed")


## get diff --------------------------------------------------------------------------------------------------------------------------------------
setwd("/media/yh/Nano2/YL/m6a/diffpeak")
re <- diffpeaksummit(peakwt, peakmut, wt_ip, wt_input, tr_ip, tr_input, recenter = F,
                     countpeak = T, countgene = F, summit = 25, extend = 50, 
                     fc_thr = 0.5, q_thr = 0.05, normalize = "gene_total")

sum(re$mut_enhance)
sum(re$WT_enhance)
sum(re$p < 0.05)/nrow(re)
sum(re$fdr < 0.05)/nrow(re)
length(unique(re[re$mut_enhance | re$WT_enhance,]$name))

source("/media/yh/Nano2/YL/m6a/diffpeak/plot.R")





