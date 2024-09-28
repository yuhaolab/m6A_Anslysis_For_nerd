setwd("/media/yh/Nano2/YL/m6a/bamnpt")

bamfiles <- list.files(path = getwd(), pattern = "[.]bam$", full.names = T)
inputbam <- bamfiles[grep("input", bamfiles)]
ipbam <- bamfiles[grep("MeRIP", bamfiles)]

ipbam <- ipbam[c(4:6, 1:3)]
inputbam <- inputbam[c(4:6, 1:3)]
id <- sub(".mark.*", "", basename(ipbam))

setwd("/media/yh/Nano2/YL/m6a/peakcalling")
source("/media/yh/Nano2/YL/m6a/peakcallingmodule/modulesFun.R")
library(GenomicFeatures)
library(rtracklayer)
library(Rsubread)
library(exomePeak)

txdb <- makeTxDbFromGFF("/media/yh/Nano2/refnopt/Arabidopsis_thaliana.TAIR10.50.gtf")
genegr <- genes(txdb)


## reads count 
setwd("/media/yh/Nano2/YL/m6a/peakcalling/peak50b25s")
peakpoint <- .makePeaksummits(genegr, extend = 50, steps = 25)

.getpeakreads(peakbed = peakpoint, ipbam = ipbam, inputbam = inputbam,
              gtf = "/media/yh/Nano2/refnopt/Arabidopsis_thaliana.TAIR10.50.gtf",
              countgene = T, strandSpecific = 0, ignoreDup = F, countMultiMappingReads = F)

