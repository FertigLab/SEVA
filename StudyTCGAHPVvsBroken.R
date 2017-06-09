library('Homo.sapiens')
library('org.Hs.eg.db')
library('GenomicRanges')
library("GSReg")
library(EBSeq)
library(limma)
library('gplots')
library(ggplot2)
library('ROCR')
library(Matrix)
library(dplyr)





HPVStatus <- read.delim("../Data/TCGA/hpv_status_tumor_site_1_14_14.txt",
                        sep = "\t",stringsAsFactors = F)
HPVStatus$Barcode <- HPVStatus$Barcode %>% gsub(replacement = "-",pattern="[.]")

HPVStatusVec <- setNames(HPVStatus$New.HPV.Status,
                          paste(HPVStatus$Barcode,"01",sep = "-"))
load("../Data/TCGA/TCGA_RSEM_processed_091615.RDa")
colnames(TCGA.RSEM) <- TCGA.RSEM %>% 
  colnames() %>% gsub(pattern="[.]",replacement="-") %>%
  substr(start = 1,stop = 15)
TCGA.RSEM.log2 <- log2(TCGA.RSEM+1)
Gene2Study.HPVNegBroken <- names(which(TCGA.RSEM.log2[,names(HPVNegBrokenVec)] %>% rowMeans()>3))






HPVNegBroken <-  read.csv("../Data/TCGA/AlteredSpliceMachineryHPVNeg.csv",
                          header = F,stringsAsFactors = F)
HPVNegBrokenVec <- setNames(HPVNegBroken$V2,HPVNegBroken$V1)
Gene2Study.HPVStatus <- names(which(TCGA.RSEM.log2[,names(HPVStatusVec)] %>% rowMeans()>3))


juncCount <- read.delim(
  "../Data/TCGA/gdac/HNSC.txt",
  header = T,sep = "\t",stringsAsFactors = F)


juncCountUnique <- juncCount[!(juncCount$Hybridization.REF %>% duplicated()),]
juncCountUnique <- juncCountUnique[-1,]

juncNames <- sapply(strsplit(x=juncCountUnique$Hybridization.REF,split = "[+]|[:]|[,]|[-]"),
       FUN = function(x) paste(x[1],":",x[2],"-",x[6],sep = ""))

rownames(juncCountUnique) <- juncNames


junc.Count.Corrected <- juncCountUnique %>% dplyr::select(contains("01A"),contains("01B"))
colnames(junc.Count.Corrected) <- junc.Count.Corrected %>% 
                                    colnames() %>%
                                    gsub(pattern="[.]",replacement = "-") %>% 
                                    substr(start = 1,stop = 15)


print("intersection of BROKEN and our TCGA:")
cat(intersect(colnames(junc.Count.Corrected),names(HPVNegBrokenVec)) %>% length())
junc.Count.Corrected.SpliceBROKEN <- junc.Count.Corrected[,names(HPVNegBrokenVec)] 
maxJuncs <- junc.Count.Corrected.SpliceBROKEN %>% apply(MARGIN = 1,FUN = max)
junc.Count.Corrected.SpliceBROKEN <- 
  junc.Count.Corrected.SpliceBROKEN[names(which(maxJuncs>5)),]

junc.Count.Corrected.SpliceBROKEN <- data.matrix(junc.Count.Corrected.SpliceBROKEN)

print("intersection of HPV's and our TCGA:")
cat(intersect(colnames(junc.Count.Corrected),names(HPVStatusVec)) %>% length())
junc.Count.Corrected.HPV <-  junc.Count.Corrected[,names(HPVStatusVec)]

junc.Count.Corrected.HPV <- data.matrix(junc.Count.Corrected.HPV)

  
save(list = c("HPVNegBrokenVec",
             "junc.Count.Corrected.SpliceBROKEN",
             "Gene2Study.HPVNegBroken"),file = "../Data/TCGA/gdac/BrokenSplice.rda")

save(list = c("HPVStatusVec",
              "junc.Count.Corrected.HPV",
              "Gene2Study.HPVStatus"),file = "../Data/TCGA/gdac/HPVPosvsHPVNeg.rda")
