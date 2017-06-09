
rename <- function(dat, oldnames, newnames) {
  datnames <- colnames(dat)
  datnames[which(datnames %in% oldnames)] <- newnames
  colnames(dat) <- datnames
  dat
}


Directory <- 
  PerturbedCount <- c(10,15,20,25)
Thresh <- 0.01
rMatsDSGenes <- vector(mode = "list",length = length(PerturbedCount))


GenesStatusrMats <- data.frame(Gene=
                                 unique(c(names(DEDS),names(DEnonDS),
                                   names(nonDEDS),names(neutralgenes))),
                               stringsAsFactors = F) %>%
                mutate(DE = Gene %in% c(names(DEDS),names(DEnonDS)),
                     DS = Gene %in% c(names(DEDS),names(nonDEDS)))



for(  i in seq_along(PerturbedCount)){
  DirectoryCur <- paste("../../../Liliana/rMATS.cuffmerge.Sim",PerturbedCount[i],
                        "/rMATS.cuffmerge.Sim",PerturbedCount[i],
                        "/MATS_output/",sep = "")
  
  Files2Read <- list.files(path=DirectoryCur,
                           pattern = "(MATS.ReadsOnTargetAndJunctionCounts.txt)$")
  
  CurDSAll <- c()
  for( j in seq_along(Files2Read)){
    SimCur <- read.delim(file = paste(DirectoryCur,Files2Read[j],sep = ""),sep = "\t",stringsAsFactors = F)
    
    
    
    SimCurDS <- SimCur %>% 
      mutate(DSStatus=PValue<Thresh & FDR<Thresh) %>% 
      group_by(geneSymbol) %>% 
      summarise(DS = any(DSStatus) ) %>% 
      as.data.frame() %>% filter(DS==TRUE) %>% 
      dplyr::select(geneSymbol) %>%
      unlist()
    
    CurDSAll <- c(CurDSAll,SimCurDS) %>% unique()
  }
  
  rMatsDSGenes[[i]] <- intersect(isos2genesvect[CurDSAll] %>% 
                                   gsub(pattern="([|].*)$",replacement="") %>%
                                   unique() ,genesChr1)
  
  GenesStatusrMats <- GenesStatusrMats %>% 
                      mutate(rMats = Gene %in% rMatsDSGenes[[i]]) %>%
                      rename("rMats",paste("rMats",PerturbedCount[[i]],sep = "")) 
  
}

sens <- GenesStatusrMats %>% 
        summarise_at(.cols = paste("rMats",PerturbedCount,sep = ""),
                     .funs = funs(sum(.==T&DS==T)/sum(DS==T)))
prec <- GenesStatusrMats %>% 
        summarise_at(.cols = paste("rMats",PerturbedCount,sep = ""),
               .funs = funs(sum(.==T&DS==T)/sum(.==T)))

sensDE <- GenesStatusrMats %>% 
  filter(DE==1) %>% 
  summarise_at(.cols = paste("rMats",PerturbedCount,sep = ""),
               .funs = funs(sum(.==T&DS==T)/sum(DS==T)))
precDE <- GenesStatusrMats %>% 
  filter(DE==1) %>% 
  summarise_at(.cols = paste("rMats",PerturbedCount,sep = ""),
               .funs = funs(sum(.==T&DS==T)/sum(.==T)))



rMatsResults <- rbind(sensitivity=sens,
                      precision=prec,
                      senesitivityDE=sensDE,
                      precisionDE = precDE)



load("../PaperSuppl/Cache/All.rda")

recallall$rMats <- sens %>% unlist()
recallDEall$rMats <- sensDE %>% unlist()
precisionall$rMats <- prec%>% unlist()
precisionDEall$rMats <- precDE%>% unlist()




# plot(x= sampNum, y= precisionall[["SEVA"]], xlab="", ylab= "", main = "Precision", type="l",col="dark red", ylim = c(0,1),pch=1)
# lines(x= sampNum, y= precisionall[["EBSEQ"]],col="blue",lty=1,pch=1)
# lines(x= sampNum, y= precisionall[["DiffSplice"]],col="green",lty=1,pch=1)
# lines(x= sampNum, y= precisionall[["rMats"]],col="grey",lty=1,pch=1)
# 
# 
# 
# lines(x= sampNum, y= precisionDEall[["SEVA"]],col="dark red",lty=2,pch=2)
# lines(x= sampNum, y= precisionDEall[["EBSEQ"]],col="blue",lty=2,pch=2)
# lines(x= sampNum, y= precisionDEall[["DiffSplice"]],col="green",lty=2,pch=2)
# lines(x= sampNum, y= precisionDEall[["rMats"]],col="grey",lty=2,pch=2)
# 
# 
# legend("bottomright",
#        legend = c("SEVA","EBSEQ","DiffSplice","rMats",
#                   "SEVA DE","EBSEQ DE","DiffSplice DE","rMats DE"),
#        col=c("dark red","blue","green","grey","dark red","blue","green","grey"), 
#        text.col = c("dark red","blue","green","grey","dark red","blue","green","grey"),
#        lty = c(1,1,1,1,2,2,2,2))
# 
# 
# plot(x= sampNum, y= recallall[["SEVA"]], xlab="", ylab= "", main = "Recall", type="l",col="dark red", ylim = c(0,1),pch=1)
# lines(x= sampNum, y= recallall[["EBSEQ"]],col="blue",lty=1,pch=1)
# lines(x= sampNum, y= recallall[["DiffSplice"]],col="green",lty=1,pch=1)
# lines(x= sampNum, y= recallall[["rMats"]],col="grey",lty=1,pch=1)
# 
# 
# lines(x= sampNum, y= recallDEall[["SEVA"]],col="dark red",lty=2,pch=2)
# lines(x= sampNum, y= recallDEall[["EBSEQ"]],col="blue",lty=2,pch=2)
# lines(x= sampNum, y= recallDEall[["DiffSplice"]],col="green",lty=2,pch=2)
# lines(x= sampNum, y= recallDEall[["rMats"]],col="grey",lty=2,pch=2)
# 
# 
# legend("bottomright",legend = c("SEVA","EBSEQ","DiffSplice","rMats","SEVA DE","EBSEQ DE","DiffSplice DE"),
#        col=c("dark red","blue","green","grey","dark red","blue","green","grey"), 
#        text.col = c("dark red","blue","green","grey","dark red","blue","green","grey"),
#        lty = c(1,1,1,1,2,2,2,2)     )

myprec <- rbind( 
  precisionall %>% 
  sapply(FUN = function(x) setNames(x,PerturbedCount)) %>% 
  melt() %>% 
  rename(c("Var1","Var2","value"),
         c("Perturbed Sample Count","Method","precision")) %>%
  mutate(Status="Differentially Spliced"),
  precisionDEall %>% 
    sapply(FUN = function(x) setNames(x,PerturbedCount)) %>% 
    melt() %>% 
    rename(c("Var1","Var2","value"),
           c("Perturbed Sample Count","Method","precision")) %>%
    mutate(Status="Differentially Spliced & Differentially Expressed"))
  

myrecall <- rbind( 
  recallall %>% 
    sapply(FUN = function(x) setNames(x,PerturbedCount)) %>% 
    melt() %>% 
    rename(c("Var1","Var2","value"),
           c("Perturbed Sample Count","Method","recall")) %>%
    mutate(Status="Differentially Spliced"),
  recallDEall %>% 
    sapply(FUN = function(x) setNames(x,PerturbedCount)) %>% 
    melt() %>% 
    rename(c("Var1","Var2","value"),
           c("Perturbed Sample Count","Method","recall")) %>%
    mutate(Status="Differentially Spliced & Differentially Expressed"))


myrecallprec <- full_join(myrecall,
                          myprec,
                          by=c("Perturbed Sample Count","Method","Status"))



pdf("../../Results/Simulation/Figure2withRMATS.pdf")
g1 <- ggplot(myrecallprec %>% 
               mutate(`Number of samples with event` =as.numeric(`Perturbed Sample Count`)),
             aes(x=`Number of samples with event`,
                 y=recall,color=Method,shape=Method,linetype=Status))+
  geom_line()+
  geom_point()+
  coord_cartesian(ylim=c(0,1))+
guides(color=FALSE,shape=FALSE,linetype=guide_legend(nrow=2))+
  theme(legend.position="bottom",
        legend.title = element_blank())


#theme(legend.position="none")

g2 <- ggplot(myrecallprec %>% 
               mutate(`Number of samples with event` =as.numeric(`Perturbed Sample Count`)),
             aes(x=`Number of samples with event`,
                 y=precision,color=Method,shape=Method,linetype=Status))+
  geom_line()+
  geom_point()+
  coord_cartesian(ylim=c(0,1))+
  guides(linetype=FALSE,
         shape=guide_legend(nrow=2),
         color=guide_legend(nrow=2))+
  theme(legend.position="bottom",
        legend.title = element_blank())#+    

print(g1)
print(g2)

print(grid.arrange(g2, g1, ncol=2))

dev.off()

pdf("../../Results/Simulation/Figure2withRMATSOnly.pdf")
print(grid.arrange(g2, g1, ncol=2))

dev.off()