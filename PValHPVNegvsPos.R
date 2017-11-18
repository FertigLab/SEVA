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


load("../Results/Real Data/TCGA/junctionPValueHPVNegPosAug.rda")
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

plotVariation(ENormal = sapply(junctionPValueHPVNegPosAug,
                               FUN = function(x) x$E1),
              ETumor = sapply(junctionPValueHPVNegPosAug,FUN = function(x) x$E2),
              DSgenes = genesHPVNegvsHPVPos,
              mainname = "DS Genes in Variation in Cancer Patients (FDR <0.01)",
              ENormalLab = paste("Variation in HPV+ (#",(HPVPos>HPVNeg) %>% sum(),">)"),
              ETumorLab = paste("Variation in HPV- (#",(HPVNeg>HPVPos) %>% sum(),">)"))
dev.off()

