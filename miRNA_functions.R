## miRNA functions
##density plot: input - cpm table, sample names,#ofsamples equals colnumber of cpm table
makeDensityPlot<-function(cpm,samples){
  require("reshape2")
  require("ggplot2")
  data=t(cpm)
  data=data.frame(sample,data)
  mdata <- melt(data, id="sample") 
  colnames(mdata) <- c("sample","miRNA","cpm")
  mdata$cpm <- log2(mdata$cpm+1)
  ggplot(mdata, aes(cpm, colour = sample)) +geom_density()+geom_vline(xintercept = log2(11),size = 1, colour = "#FF3721",
                                                                      linetype = "dashed")+labs(x="log2(CPM+1)")#+geom_text(x=4.5,y=0.5,color="black",label="CPM=10")
  ggsave(file="densityPlot.eps")
}