ex.test <- function(dgelist,mean=NULL,g1,g2,fc=1.5,p=0.05,replicates=T,disp=0.16,pcounts){
  #View(mean)
  if(replicates==T)
    {et <- exactTest(dgelist,pair=c(g1,g2))}
  else
    {print("no replicates in given data")
    et <- exactTest(dgelist,pair=c(g1,g2),dispersion=disp)
    }
  
  list <-topTags(et,n=nrow(dgelist))
  dim(list)
  View(list)
  if(length(mean)!=0){
    print("merge exacttest result and counts table")
    fls <- merge(list,mean,by="row.names")
    row.names(fls) <- fls$Row.names
    fls <- fls[,-1]
  }
 
  else{
    fls <- list$table
  }
 
  
      
  print(nrow(fls))
  print("is number of rows of fls")
  FC <-abs(fls$logFC) >= log2(fc)
  #print(summary(FC))
  
  #default.set = FC & fls$PValue <= p
  #default.set = FC & fls$FDR <= p
  default.set= rep("no change",nrow(fls))
  #up <- fls$logFC>=log2(fc) & fls$FDR<=p
  up <- fls$logFC>=log2(fc) & fls$PValue<=p
  #down <- fls$logFC<=-log2(fc) & fls$FDR<=p
  down <- fls$logFC<=-log2(fc) & fls$PValue<=p
  default.set[up] <- "up-regulated"
  default.set[down] <- "down-regulated"
  
  
  
  print(summary(default.set))
  #cat(length(default.set))
  #View(fls[order(row.names(fls)),])
  fls1 <- data.frame(fls, default.set)
  #View(fls1[order(row.names(fls1)),])
#   pattern1<-paste("_",g1,"_|","_",g2,"_",sep="")
#   print(pattern1)
#   foo<-c(group,group)
#   keep <- grep(pattern1,foo)
#  print(keep)
#  pcounts2<-pcounts[,keep]
  
#  View(pcounts2)
  fls2 <- merge(fls1,pcounts,by="row.names")
  #View(fls2[order(row.names(fls2)),])
  #View(fls2)
  #print (paste("There are", sum(fls2$default.set & fls2$logFC>0), "genes up regulated and", sum(fls2$default.set & fls2$logFC<=0), "genes down regulated in Group", g2,  "vs Group", g1,"."))
  print (paste("There are", sum(fls2$default.set=="up-regulated"), "genes up regulated and", sum(fls2$default.set =="down-regulated"), "genes down regulated in Group", g2,  "vs Group", g1,"."))

  print (paste("There are", length(default.set[default.set!="no change"]), "genes meeting the default criteria which are differntially expressed in Group", g2,  "vs Group", g1,"."))
  colnames(fls2)[1]<-"Gene_Symbol"
  fls2<-fls2[,-3]   ### remove column of logCPM
  write.table(fls2,paste(g2,"_vs_",g1,".txt",sep=""),sep="\t",row.names=F,quote=F)
  print(paste("--------Done analysis for",g2,"vs",g1,"---------------"))
  fls2
}