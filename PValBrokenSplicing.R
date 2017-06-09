load("../Data/TCGA/gdac/BrokenSplice.rda")
#We divide the TCGA data to batches of 1000 genes to make it manageable
junctionPValueBROKENAug <- c()
gene_partition <- c(seq(1,length(Gene2Study.HPVNegBroken),1000),
                    length(Gene2Study.HPVNegBroken))

for( i in seq_along(gene_partition[-1])){
  
  junctionPValueBROKEN <- GSReg.SEVA(juncExprs = junc.Count.Corrected.SpliceBROKEN,  
                                     phenoVect=as.factor(HPVNegBrokenVec),
                                   sparse = T, 
                                   verbose = F,
                                   GenestoStudy = 
                                   Gene2Study.HPVNegBroken[gene_partition[i]:gene_partition[i+1]])
  
  
  print(i)
  gc()
  junctionPValueBROKENAug <- c(junctionPValueBROKENAug,junctionPValueBROKEN)
  save(file = "../Results/Real Data/TCGA/junctionPValueBROKENAug.rda",
       list = "junctionPValueBROKENAug")
}


pdf("../Results/Real Data/TCGA/BrokenVsUnBrokenBF.pdf")
genesHPVNegvsHPVPos <- 
  names(which(sapply(junctionPValueBROKENAug,function(x) x$pvalueTotal)<
                0.01/length(junctionPValueBROKENAug)))

UnBroken <- sapply(junctionPValueBROKENAug[genesHPVNegvsHPVPos],
             FUN = function(x) x$E1)
Broken <- sapply(junctionPValueBROKENAug[genesHPVNegvsHPVPos],
             FUN = function(x) x$E2)


plot(x=UnBroken,y=Broken,
     xlim=c(0,0.31),
     ylim=c(0,0.31),
     main="DS Genes in Variation in HPV- Patients (Bonferroni <0.01)",
     xlab = paste("Variation in Unaltered Spliced Machinery (#",(UnBroken>Broken) %>% sum(),">)"),
     ylab = paste("Variation in Altered Spliced Machinery (#",(UnBroken<Broken) %>% sum(),">)"))

lines(x=c(-.05,0.4),y=c(-0.05,0.4),
      type = "l",
      lty = 2,
      col="dark red")

dev.off()

pdf("../Results/Real Data/TCGA/BrokenVsUnBrokenFDR.pdf")


genesHPVNegvsHPVPos <- 
  (p.adjust(sapply(junctionPValueBROKENAug,
          function(x) x$pvalueTotal),method = "BH")<0.01) %>% 
                          which() %>% names()

UnBroken <- sapply(junctionPValueBROKENAug[genesHPVNegvsHPVPos],
                   FUN = function(x) x$E1)
Broken <- sapply(junctionPValueBROKENAug[genesHPVNegvsHPVPos],
                 FUN = function(x) x$E2)

plot(x=UnBroken,y=Broken,
     xlim=c(0,0.31),
     ylim=c(0,0.31),
     main="DS Genes in Variation in HPV- Patients (FDR <0.01)",
     xlab = paste("Variation in Unaltered Spliced Machinery (#",(UnBroken>Broken) %>% sum(),">)"),
     ylab = paste("Variation in Altered Spliced Machinery (#",(UnBroken<Broken) %>% sum(),">)"))

lines(x=c(-.05,0.4),y=c(-0.05,0.4),
      type = "l",
      lty = 2,
      col="dark red")
dev.off()

