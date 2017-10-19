exon_reads_count <- function(file,exon_anno,read_length=51){
  require("GenomicAlignments")
  require("GenomicRanges")
  cords <- c(exon_anno$start,exon_anno$end)
  start<-min(cords);end<-max(cords)
  chr <- as.character(exon_anno$chr[1])
  data <- readGAlignments(file, use.names=TRUE, param = ScanBamParam(which=GRanges(chr, IRanges(start, end))))
  data_cov <- coverage(data)
  ann_chr <- exon_anno
  #ann_chr$length <- ann_chr$end - ann_chr$start+1  
  ann_chr$cov <- viewSums(Views(data_cov[[chr]], start=ann_chr$start, end=ann_chr$end))
  ann_chr$count <- round(viewSums(Views(data_cov[[chr]], start=ann_chr$start, end=ann_chr$end))/read_length,1)  ## 100 is the read length
  exon_count <- ann_chr
  exon_count_sum <- sum(exon_count$count)
  exon_count_sum
}