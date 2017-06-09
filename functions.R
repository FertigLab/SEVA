#### functions for SEVA paper
#### Bahman Afsari


library('Homo.sapiens')
library('org.Hs.eg.db')
library('GenomicRanges')
library("GSReg")


### DEgenes
findDEGenes <- function(geneexp,phenoVect,pvaluecorrected=0.01){
  geneexpr <- as.matrix(geneexp)

  Sizes = MedianNorm(geneexpr)
  EBOut = EBTest(Data = geneexpr,
               Conditions = as.factor(phenoVect),
               sizeFactors = Sizes, maxround = 5)


  PP = GetPPMat(EBOut)
  DEGenes_EBesq <- as.character(names(which(PP[,"PPDE"]>1-pvaluecorrected)))
  return(list(DEGenes_EBesq=DEGenes_EBesq,pvaluescorrected=PP[,"PPDE"]))
}


### EBSEQ isoform DS
findDSEBSEQ<-function(isos,phenoVect, isos2genesvectsimplified, pvaluecorrected= 0.05 )
{
  isoexpressed <- names(which(apply(isos,MARGIN = 1, FUN = sum)>0.0000001))
  Sizes = MedianNorm(isos[isoexpressed,])
  IsoEBOut = EBTest(Data = isos[isoexpressed,],
                    Conditions = as.factor(phenoVect),
                    sizeFactors = Sizes, maxround = 5)
  
  
  PPISO = GetPPMat(IsoEBOut)
  ebseqpvalueGenes <- tapply(  X = PPISO[,"PPDE"] ,
                               INDEX = isos2genesvectsimplified[rownames(PPISO)],
                               FUN = function(x) mean(x,na.rm=T))
  DSEBseq <- names(which(ebseqpvalueGenes>1-pvaluecorrected))
  return(list(DSEBseq=DSEBseq, pvaluescorrected=ebseqpvalueGenes))
  
}

### SEVA
### Make the matrix of overlaps
OverlapJunctionMatrices<-function(junc.RPM, GenestoStudy=NULL){
  ##loading genes
  gnAll <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  gSymbol <- select(org.Hs.eg.db,keys=as.character(gnAll$gene_id),
                    columns=c('SYMBOL'),keytype='ENTREZID')
  gnAll$SYMBOL <- gSymbol$SYMBOL
  
  if(is.null(GenestoStudy)){
    gn <- gnAll
  }else{
    GenestoStudyId <- which(is.element(gnAll$SYMBOL,GenestoStudy))
    gn <- gnAll[GenestoStudyId]
  }
  
  
  ### finding junction locations and put them in the size
  z <- strsplit(sub('-',':',rownames(junc.RPM)),':')
  mychr <- sapply(X=z,FUN = function(x) x[1])
  mystart <- sapply(X=z,FUN = function(x) x[2])
  myend <- sapply(X=z,FUN = function(x) x[3])
  
  ### puting junction in GRanges format
  juncRanges <- as.data.frame(strsplit(sub('-',':',rownames(junc.RPM)),':'),stringsAsFactors = FALSE)
  junctionsGRanges <- GRanges(seqnames = Rle(mychr), 
                              ranges = IRanges(start = as.numeric(mystart), end = as.numeric(myend)))
  
  ### finding junction for each gene
  overlapJunction <- findOverlaps(junctionsGRanges,gn)
  
  #making overlap matrix
  juncnames <- rownames(junc.RPM)
  genesJunction <- tapply(rownames(junc.RPM)[queryHits(overlapJunction)],
                          gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)
  
  genesJunctionInd <- tapply(queryHits(overlapJunction),
                             gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)
  MyRest <- sapply(X = genesJunctionInd, function(x) { 
    y <- junctionsGRanges[x];
    w <- findOverlaps(y,y);
    mymat <- Matrix(0,nrow = length(y), ncol = length(y),
                    dimnames = list(juncnames[x],juncnames[x]),sparse = T);
    mymat[queryHits(w)+ length(x)*(subjectHits(w)-1)] <- 1
    return(mymat)
  })
  return(list(Rest = MyRest,genesJunction = genesJunction))
}
### Main function
SEVA.meangeneFilter <- function(junc.RPM,phenoVect,
                               geneexpr,minmeanloggeneexp= 3, GenestoStudy = NULL){
  
  z <- OverlapJunctionMatrices(junc.RPM, GenestoStudy)
  MyRest <- z$Res 
  genesJunction <- z$genesJunction
  geneexpressed <- names(which(apply(log2(geneexpr+1),MARGIN = 1,FUN = mean)>minmeanloggeneexp))
  
  junctionPValue <- GSReg.GeneSets.EVA(geneexpres = junc.RPM[,names(phenoVect)],phenotypes = as.factor(phenoVect), 
                                       minGeneNum = 2,
                                       pathways = genesJunction[intersect(names(which(sapply(genesJunction ,length) >2)),
                                                                          geneexpressed)],
                                       distFunc = GSReg.kendall.tau.distance.Restricted.Sparse,
                                       distparamPathways = MyRest  )
  return(junctionPValue)
}

### You can call this function to filter the junctions that whose log2(expression+1) is 
# less than alpha* log2(its gene expression + 1). The function returns the junction names
SEVA.meangenemaxjunc.Filter <- function(junc.RPM, geneexpr, alpha =0.1, GenestoStudy=NULL){
  ##loading genes
  gnAll <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  gSymbol <- select(org.Hs.eg.db,keys=as.character(gnAll$gene_id),
                    columns=c('SYMBOL'),keytype='ENTREZID')
  gnAll$SYMBOL <- gSymbol$SYMBOL
  
  if(is.null(GenestoStudy)){
    gn <- gnAll
  }else{
    GenestoStudyId <- which(is.element(gnAll$SYMBOL,GenestoStudy))
    gn <- gnAll[GenestoStudyId]
  }
  
  
  
  ### finding junction locations and put them in the size
  z <- strsplit(sub('-',':',rownames(junc.RPM)),':')
  mychr <- sapply(X=z,FUN = function(x) x[1])
  mystart <- sapply(X=z,FUN = function(x) x[2])
  myend <- sapply(X=z,FUN = function(x) x[3])
  
  ### puting junction in GRanges format
  juncRanges <- as.data.frame(strsplit(sub('-',':',rownames(junc.RPM)),':'),stringsAsFactors = FALSE)
  junctionsGRanges <- GRanges(seqnames = Rle(mychr), 
                              ranges = IRanges(start = as.numeric(mystart), end = as.numeric(myend)))
  
  ### finding junction for each gene
  overlapJunction <- findOverlaps(junctionsGRanges,gn)
  
  
  
  #making overlap matrix
  juncnames <- rownames(junc.RPM)
  #   genesJunction <- tapply(rownames(junc.RPM)[queryHits(overlapJunction)],
  #                           gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)
  #   
  #   genesJunctionInd <- tapply(queryHits(overlapJunction),
  #                              gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)
  #   
  
  #genes junction with their index in junc.RPM
  genesJunction <- tapply(gn$SYMBOL[subjectHits(overlapJunction)],
                          juncnames[queryHits(overlapJunction)],
                          function(x) x[1])
  
  
  
  meanloggeneexpr <- apply(log2(geneexpr+1),MARGIN = 1,FUN = mean)
  juncmax <-  apply(log2(junc.RPM+1),MARGIN = 1, max)
  
  #geneexpressed <- names(which(meanloggeneexpr>minmeanloggeneexp))
  juncexpressed <- names(which(juncmax[names(genesJunction)]>0.1*meanloggeneexpr[genesJunction]))
  
  #geneexpressed <- names(which(apply(log2(geneexpr+1),MARGIN = 1,FUN = mean)>minmeanloggeneexp))
  
  return(juncexpressed)
}



#junctions to filter
# SEVA.meangenemaxjunc.Filter <- function(junc.RPM,phenoVect,
#                                 geneexpr,minmeanloggeneexp= 3, alpha =0.1, GenestoStudy = NULL){
#   ##loading genes
#   gnAll <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#   gSymbol <- select(org.Hs.eg.db,keys=as.character(gnAll$gene_id),
#                     columns=c('SYMBOL'),keytype='ENTREZID')
#   gnAll$SYMBOL <- gSymbol$SYMBOL
#   
#   if(is.null(GenestoStudy)){
#     gn <- gnAll
#   }else{
#     GenestoStudyId <- which(is.element(gnAll$SYMBOL,GenestoStudy))
#     gn <- gnAll[GenestoStudyId]
#   }
#   
#   
#   
#   
#   ### finding junction locations and put them in the size
#   z <- strsplit(sub('-',':',rownames(junc.RPM)),':')
#   mychr <- sapply(X=z,FUN = function(x) x[1])
#   mystart <- sapply(X=z,FUN = function(x) x[2])
#   myend <- sapply(X=z,FUN = function(x) x[3])
#   
#   ### puting junction in GRanges format
#   juncRanges <- as.data.frame(strsplit(sub('-',':',rownames(junc.RPM)),':'),stringsAsFactors = FALSE)
#   junctionsGRanges <- GRanges(seqnames = Rle(mychr), 
#                               ranges = IRanges(start = as.numeric(mystart), end = as.numeric(myend)))
#   
#   ### finding junction for each gene
#   overlapJunction <- findOverlaps(junctionsGRanges,gn)
#   
#   
#   
#   #making overlap matrix
#   juncnames <- rownames(junc.RPM)
# #   genesJunction <- tapply(rownames(junc.RPM)[queryHits(overlapJunction)],
# #                           gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)
# #   
# #   genesJunctionInd <- tapply(queryHits(overlapJunction),
# #                              gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)
# #   
#   
#   #genes junction with their index in junc.RPM
#   genesJunction <- tapply(queryHits(overlapJunction),
#                              gn$SYMBOL[subjectHits(overlapJunction)],
#                              function(x) 
#                               {
#                                   names(x)<-rownames(junc.RPM)[x];
#                                   return(x);
#                               })
#   
#   
#   meanloggeneexpr <- apply(log2(geneexpr+1),MARGIN = 1,FUN = mean)
#   juncmax <-  apply(log2(junc.RPM+1),MARGIN = 1, max)
#   
#   geneexpressed <- names(which(meanloggeneexpr>minmeanloggeneexp))
#   
#   
#   genesJunctionExpressed <- sapply(intersect(intersect(names(genesJunction),
#                                                        names(meanloggeneexpr)),geneexpressed), 
#          FUN =  function(x) genesJunction[[x]][which(juncmax[genesJunction[[x]]]>meanloggeneexpr[x]*0.1)])
#   
#   
#   FinalGenes <- names(which(sapply(genesJunctionExpressed,length)>2))
#   
#   MyRest <- sapply(X = genesJunctionExpressed[FinalGenes], function(x) { 
#     y <- junctionsGRanges[x];
#     w <- findOverlaps(y,y);
#     mymat <- matrix(0,nrow = length(y), ncol = length(y),dimnames = list(juncnames[x],juncnames[x]));
#     mymat[queryHits(w)+ length(x)*(subjectHits(w)-1)] <- 1
#     return(mymat)
#   })
#   #geneexpressed <- names(which(apply(log2(geneexpr+1),MARGIN = 1,FUN = mean)>minmeanloggeneexp))
#   
#   junctionPValue <- GSReg.GeneSets.EVA(geneexpres = junc.RPM[,names(phenoVect)],phenotypes = as.factor(phenoVect), 
#                                        minGeneNum = 2,
#                                        pathways = sapply(genesJunction[FinalGenes],names),
#                                        distFunc = GSReg.kendall.tau.distance.Restricted,
#                                        distparamPathways = MyRest  )
#   return(junctionPValue)
# }
# 


### plot variation diagram and if VennMatrix is available it plots venn diagram as well
#Venn columns must be DE, EBSEQ, DiffSplice and SEVA
plotVariation <- function(ENormal,ETumor, DSgenes,mainname = deparse(substitute(DSgenes)), VennMatrix){
  #mainname with number of identified genes with higher variation in tumore and in cancer
  if(missing(VennMatrix)){
    AllGenes <- names(ENormal)
  }else{
    AllGenes <- intersect(names(ENormal),rownames(VennMatrix))
  }
  
  DSgenesIntersect <- intersect(AllGenes,DSgenes)
  
  mainename_variation_count <- paste0(mainname, " # Tumor>[Normal>]", #main naime
                                      sum(ETumor[DSgenesIntersect] > ENormal[DSgenesIntersect]), "[",# gene number with higher variation in tumor
                                      sum(ETumor[DSgenesIntersect] < ENormal[DSgenesIntersect]), "]") # gene number with higher variation in cancer
  
  
  
  
  Erange <- range(c(ENormal,ETumor)) #range of variation
  plot(x=ENormal[setdiff(names(ENormal),DSgenesIntersect)], #plot non-DS
       y=ETumor[setdiff(names(ETumor),DSgenesIntersect)], main = mainename_variation_count,
       col="light blue", xlab = "Normal Variation", ylab="Tumor Variation",
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

GSReg.kendall.tau.distance.Restricted.Sparse <- function(V, RestMat){
  #Checking if V is a numeric matrix 
  
  GSReg:::GSReg.Check.input(V)
  
  inter = intersect(rownames(V),rownames(RestMat))
  if( length(inter)< 2)
    stop("No common restrictions.")
  V = V[inter,]
  RestMat <- as(RestMat[inter,inter],"matrix")
  
  n <- dim(V)[1]
  m <- dim(V)[2]
  
  d <- .C("kendalltaudistRestricted",
          as.double(V),
          as.integer(dim(V)[1]),
          as.integer(dim(V)[2]),
          as.integer(RestMat),
          as.double(matrix(data= 0 ,ncol=m,nrow=m)))
  dist <- d[[5]]
  dim(dist) <- c(m,m)
  rownames(dist) <- colnames(V)
  colnames(dist) <- colnames(V)
  return(dist/(sum(RestMat)+0.000001))
}

