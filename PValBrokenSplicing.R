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

load("../Results/Real Data/TCGA/junctionPValueBROKENAug.rda")
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


plotVariation <- function(ENormal,ETumor,
                          DSgenes,mainname = deparse(substitute(DSgenes)),
                          ENormalLab = "Normal Variation",
                          ETumorLab = "Tumor Variation",
                          VennMatrix,AddtoMain = F){
  #mainname with number of identified genes with higher variation in tumore and in cancer
  if(missing(VennMatrix)){
    AllGenes <- names(ENormal)
  }else{
    AllGenes <- intersect(names(ENormal),rownames(VennMatrix))
  }
  
  DSgenesIntersect <- intersect(AllGenes,DSgenes)
  
  if(AddtoMain){
    mainename_variation_count <- paste0(mainname, " # Tumor>[Normal>]", #main naime
                                        sum(ETumor[DSgenesIntersect] > ENormal[DSgenesIntersect]), "[",# gene number with higher variation in tumor
                                        sum(ETumor[DSgenesIntersect] < ENormal[DSgenesIntersect]), "]") # gene number with higher variation in cancer
  }else{
    mainename_variation_count <- mainname
  }
  
  
  
  Erange <- range(c(ENormal,ETumor)) #range of variation
  plot(x=ENormal[setdiff(names(ENormal),DSgenesIntersect)], #plot non-DS
       y=ETumor[setdiff(names(ETumor),DSgenesIntersect)], main = mainename_variation_count,
       col="light blue", xlab =ENormalLab , ylab=ETumorLab,
       xlim = Erange,ylim = Erange,pch = 18)
  lines(x=ENormal[DSgenesIntersect],y=ETumor[DSgenesIntersect], #plot DS
        col="dark red",type = "p",pch = 19)
  lines(x=Erange,y=Erange,col="black",type = "l", lty = 2,lwd = 2) #45 degree line
  legend("topleft", legend = c("DS", "non-DS"),pch = c(20,18), #legend
         col = c("dark red","light blue"),
         text.col = c("dark red","light blue"))
  
  if(!missing(VennMatrix)){
    VennMatrixCopy <- VennMatrix
    VennMatrixCopy[,"SEVA"] <- 0
    VennMatrixCopy[DSgenesIntersect,"SEVA"] <- 1
    vennDiagram(VennMatrixCopy[,c("SEVA","EBSEQ","DE")])
    vennDiagram(VennMatrixCopy[,c("SEVA","DiffSplice","DE")])  
  }
  
}


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

plotVariation(ENormal = sapply(junctionPValueBROKENAug,
                               FUN = function(x) x$E1),
              ETumor = sapply(junctionPValueBROKENAug,FUN = function(x) x$E2),
              DSgenes = genesHPVNegvsHPVPos,
              mainname = "DS Genes in Variation in HPV- Patients (FDR <0.01)",
              ENormalLab =  paste("Variation in Unaltered Spliced Machinery (#",(UnBroken>Broken) %>% sum(),">)"),
              ETumorLab = paste("Variation in Altered Spliced Machinery (#",(UnBroken<Broken) %>% sum(),">)"))
dev.off()

