setwd("//isi-dcnl/user_data/Seq/cwarden/HTseq_Files/")
genome = "hg19"
annotation.file = paste("TxDb_",genome,"_exon_annotations.txt",sep="")
gene.length.file = paste("TxDb_",genome,"_gene_length_and_position.txt",sep="")

annotation.table = read.table(annotation.file, sep="\t", head=T)

gene.symbol = as.character(levels(as.factor(as.character(annotation.table$symbol))))
gene.chr = annotation.table$chr[match(gene.symbol, annotation.table$symbol)]
gene.start= tapply(as.numeric(annotation.table$start), as.character(annotation.table$symbol), min)
gene.end = tapply(as.numeric(annotation.table$end), as.character(annotation.table$symbol), max)
gene.strand = annotation.table$strand[match(gene.symbol, annotation.table$symbol)]
gene.length = tapply(as.numeric(annotation.table$width), as.character(annotation.table$symbol), sum)

output.table = data.frame(symbol = gene.symbol, chr=gene.chr, start=gene.start, end=gene.end, strand=gene.strand, length=gene.length)
write.table(output.table, gene.length.file, sep="\t", row.names=F, quote=F)