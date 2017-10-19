get.symbolIDs <- function(id,id.type){
  
  # ************************************************
  
  # get.symbolIDs function programmed by 
  # Benjamin Tovar | February 25, 2014
  
  # INSTRUCTIONS
  # id = vector of the original IDs, for example: 
  # 	c("ENSG00000121410","ENSG00000175899")
  # id.type = type of the original IDs, in the example is ensembl
  
  # NOTE: only refseq, ensembl, uniprot, unigene are supported,
  # this function depends on the Bioconductor package org.Hs.eg.db
  # that can be installed:
  
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("org.Hs.eg.db")
  
  # ************************************************
  # # USAGE EXAMPlE: ENSEMBL
  # 	require(org.Hs.eg.db)
  # 	ensembl <- toTable(org.Hs.egENSEMBL)
  # 	id  <- ensembl[1:100,2]
  # 	id.type  <- "ensembl"
  # 	res <- get.symbolIDs(id,id.type)
  
  # # USAGE EXAMPlE: UNIPROT
  #   require(org.Hs.eg.db)
  # 	uniprot <- toTable(org.Hs.egUNIPROT)
  # 	id  <- uniprot[1:100,2]
  # 	id.type  <- "uniprot"
  # 	res <- get.symbolIDs(id,id.type)
  
  # # USAGE EXAMPlE: REFSEQ
  # require(org.Hs.eg.db)
  # refseq.id <- toTable(org.Hs.egREFSEQ)
  # id  <- refseq.id[1:100,2]
  # id.type  <- "refseq"
  # res <- get.symbolIDs(id,id.type)	
  
  # # USAGE EXAMPlE: UNIGENE
  # require(org.Hs.eg.db)
  # unigene <- toTable(org.Hs.egUNIGENE)
  # id  <- unigene[1:100,2]
  # id.type  <- "unigene"
  # res <- get.symbolIDs(id,id.type)	
  
  # LOAD THE ANNOTATION INFORMATION
  cat("1) Loading annotation library",date(),"\n")   	
  require(org.Hs.eg.db)
  cat("2) Loading annotation for symbol IDs",date(),"\n")   	    
  symbol <- toTable(org.Hs.egSYMBOL)
  # IF THE ORIGINALS IDS = ENSEMBL
  if(id.type=="ensembl"){
    cat("3) Loading annotation for ensembl IDs",date(),"\n")   	
    annotation <- toTable(org.Hs.egENSEMBL)
    # extract the indexes of the original database 
    index <- as.numeric(sapply(id, function(x) which(annotation[,2]==x)[1]))
    # extract the gene_ids
    gene_id.index <- as.numeric(annotation[index,1])
    # parse the indexes back to the symbol IDs 
    index <- as.numeric(sapply(gene_id.index, function(x) which(symbol[,1]==x)))
    # extract the IDs
    symbolIDs <- symbol[index,2]
    # export the output
    return(symbolIDs)
  }
  # IF THE ORIGINALS IDS = UNIPROT
  if(id.type=="uniprot"){
    cat("3) Loading annotation for uniprot IDs",date(),"\n")   
    annotation <- toTable(org.Hs.egUNIPROT)
    # extract the indexes of the original database 
    index <- as.numeric(sapply(id, function(x) which(annotation[,2]==x)[1]))
    # extract the gene_ids
    gene_id.index <- as.numeric(annotation[index,1])
    # parse the indexes back to the symbol IDs 
    index <- as.numeric(sapply(gene_id.index, function(x) which(symbol[,1]==x)))
    # extract the IDs
    symbolIDs <- symbol[index,2]
    # export the output
    return(symbolIDs)
  }
  # IF THE ORIGINALS IDS = REFSEQ
  if(id.type=="refseq"){
    cat("3) Loading annotation for refseq IDs",date(),"\n")   
    annotation <- toTable(org.Hs.egREFSEQ)
    # extract the indexes of the original database 
    index <- as.numeric(sapply(id, function(x) which(annotation[,2]==x)[1]))
    # extract the gene_ids
    gene_id.index <- as.numeric(annotation[index,1])
    # parse the indexes back to the symbol IDs 
    index <- as.numeric(sapply(gene_id.index, function(x) which(symbol[,1]==x)))
    # extract the IDs
    symbolIDs <- symbol[index,2]
    # export the output
    return(symbolIDs)
  }
  # IF THE ORIGINALS IDS = UNIGENE
  if(id.type=="unigene"){
    cat("3) Loading annotation for unigene IDs",date(),"\n")   
    annotation <- toTable(org.Hs.egUNIGENE)
    # extract the indexes of the original database 
    index <- as.numeric(sapply(id, function(x) which(annotation[,2]==x)[1]))
    # extract the gene_ids
    gene_id.index <- as.numeric(annotation[index,1])
    # parse the indexes back to the symbol IDs 
    index <- as.numeric(sapply(gene_id.index, function(x) which(symbol[,1]==x)))
    # extract the IDs
    symbolIDs <- symbol[index,2]
    # export the output
    return(symbolIDs)
  }		
  cat("** ERROR: DATABASE TYPE NOT SELECTED | TRY: ensembl OR unigene OR uniprot OR refseq **",date(),"\n") 
  return(NA)
}