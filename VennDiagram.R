library("VennDiagram", lib.loc="C:/Program Files/R/R-3.1.2/library")
venn.diagram(list(A=shCD133.up,B=shGlis3.up),filename = "overlap-up-genes.tiff",fill=c("yellow","blue"),category.names = c("shCD133","shGlis3"))
