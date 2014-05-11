outToExcel<-function(list,file="result.xls",sheetnames){
  if(length(list)!=length(sheetnames)) stop("unequal length of output elements and sheets")
  n = length(sheetnames)
  library(xlsx)
  if(file %in% list.files()) stop("output file already exists")
  for (i in 1:n){
    write.xlsx2(list[[i]],file,sheetName=sheetnames[i],append=T)
  }
}