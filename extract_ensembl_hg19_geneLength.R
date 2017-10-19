### extract total non-overlapping exon length per gene with R bioconductor
### 170515---- Lu Yang
library(GenomicFeatures)
# First, import the GTF-file that you have also used as input for htseq-count
ref.file <- "//isi-dcnl/user_data/Seq/Lu/ref/Homo_sapiens.GRCh37.87.protein_coding.and.lncRNA.gtf"
hg19.Ensembel.txdb <- makeTxDbFromGFF(ref.file, format="gtf")
# then collect the exons per gene id
#id2name(hg19.Ensembel.txdb,)
exons.list.per.gene <- exonsBy(hg19.Ensembel.txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
## convert ensembel id to gene name
require(biomaRt)
# Collect ensembl IDs from annot.table before converting to normal gene names.
annot.table <- t(data.frame(exonic.gene.sizes));colnames(annot.table) <- "gene_length"
ensembl_ids <- rownames(annot.table)
# Prepare gene table with some simple caching to avoid stressing the Ensembl server by many repeated runs
genes.table = NULL
###R code for fetching gene names with biomaRt ###
if (!file.exists("cache.genes.table")){
  message("Retrieving genes table from Ensembl...")
  mart <- useMart("ensembl")
  grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", dataset="hsapiens_gene_ensembl")
  #listDatasets(mart=mart)
  #mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
  #genes.table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description"), values= ensembl_ids, mart= mart)
  genes.table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "description"), values= ensembl_ids, mart= grch37)
  #save(genes.table, file= "cache.genes.table")
  
} else {
  load("cache.genes.table")
  message("Reading gene annotation information from cache file: cache/cache.genes.table Remove the file if you want to force retrieving data from Ensembl")
  
}
# Merging two tables,
annot.table <- merge(x = annot.table, y = genes.table, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, all.y = F )
colnames(annot.table)[1]<-"ensembl_id"
save(annot.table,file="ensembl_hg19_geneLength.rd")
write.table(annot.table,file="ensembl_hg19_geneLength.txt",sep="\t",row.names = F)


