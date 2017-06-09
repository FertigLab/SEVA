load("../Data/TCGA/gdac/HPVPosvsHPVNeg.rda")
#We divide the TCGA data to batches of 1000 genes to make it manageable
junctionPValueHPVNegPosAug <- c()
gene_partition <- c(seq(1,length(Gene2Study.HPVStatus),1000),
                    length(Gene2Study.HPVStatus))

for( i in seq_along(gene_partition[-1])){
  
  junctionPValueHPVNegPos <- GSReg.SEVA(juncExprs = junc.Count.Corrected.HPV,  
                                        phenoVect=as.factor(HPVStatusVec),
                                        sparse = T, 
                                        verbose = F,
                                        GenestoStudy = 
                                          Gene2Study.HPVStatus[gene_partition[i]:gene_partition[i+1]])
  
  
  print(i)
  gc()
  junctionPValueHPVNegPosAug <- c(junctionPValueHPVNegPosAug,junctionPValueHPVNegPos)
  save(file = "../Results/Real Data/TCGA/junctionPValueHPVNegPosAug.rda",
       list = "junctionPValueHPVNegPosAug")
}


pdf("../Results/Real Data/TCGA/HPVNegPosBF.pdf")
genesHPVNegvsHPVPos <- 
  names(which(sapply(junctionPValueHPVNegPosAug,function(x) x$pvalueTotal)<
                0.01/length(junctionPValueHPVNegPosAug)))

prin(levels(as.factor(HPVStatusVec)))

HPVNeg <- sapply(junctionPValueHPVNegPosAug[genesHPVNegvsHPVPos],
             FUN = function(x) x$E1)
HPVPos <- sapply(junctionPValueHPVNegPosAug[genesHPVNegvsHPVPos],
             FUN = function(x) x$E2)


plot(x=HPVNeg,y=HPVPos,
     xlim=c(0,0.31),
     ylim=c(0,0.31),
     main="DS Genes in Variation in Cancer Patients (Bonferroni <0.01)",
     ylab = paste("Variation in HPV- (#",(HPVNeg>HPVPos) %>% sum(),">)"),
     xlab = paste("Variation in HPV+ (#",(HPVPos>HPVNeg) %>% sum(),">)"))

lines(x=c(-.05,0.4),y=c(-0.05,0.4),
      type = "l",
      lty = 2,
      col="dark red")

dev.off()

pdf("../Results/Real Data/TCGA/HPVNegPosFDR.pdf")


genesHPVNegvsHPVPos <- 
  (p.adjust(sapply(junctionPValueHPVNegPosAug,
          function(x) x$pvalueTotal),method = "BH")<0.01) %>% 
                          which() %>% names()

HPVNeg <- sapply(junctionPValueHPVNegPosAug[genesHPVNegvsHPVPos],
                 FUN = function(x) x$E1)
HPVPos <- sapply(junctionPValueHPVNegPosAug[genesHPVNegvsHPVPos],
                 FUN = function(x) x$E2)


plot(x=HPVNeg,y=HPVPos,
     xlim=c(0,0.31),
     ylim=c(0,0.31),
     main="DS Genes in Variation in Cancer Patients (FDR <0.01)",
     ylab = paste("Variation in HPV- (#",(HPVNeg>HPVPos) %>% sum(),">)"),
     xlab = paste("Variation in HPV+ (#",(HPVPos>HPVNeg) %>% sum(),">)"))

lines(x=c(-.05,0.4),y=c(-0.05,0.4),
      type = "l",
      lty = 2,
      col="dark red")
dev.off()

