
########### richness, diversity and evenness in this function
Evenness <- function(mat) {
  require(vegan)
  H1 <- diversity(mat)
  R <- specnumber(mat)
  J <- H1 / log(R)
  hrj <- data.frame( H1, R, J )
  invisible(hrj)
}


################## makes taxonomy data.frame from taxon vector, sep='; '
makeTaxo <- function(taxo.in=tax.np, otu.table=pb.3500) {
	taxo.in <- as.character(taxo.in)
	tax.tmp.ls <- strsplit(taxo.in, split='; ')
	tax.lengths <- unlist(lapply(tax.tmp.ls, length))
	max(tax.lengths)
	tax.tmp.ls[[1]][1]
	
	# test
	# x <- c('a; b; c; d', 'e; f; g', 'i; j')
	# x2 <- strsplit(x, '; ')
	# x3 <- data.frame(one=sapply(x2, function(x){x[1]}),
					 # two=sapply(x2, function(x){x[2]}),
					 # three=sapply(x2, function(x){x[3]}),
					 # four=sapply(x2, function(x){x[4]}))
	# x3
	# x3$four <- as.character(x3$four)
	# x3$four[which(is.na(x3$four))] <- 'h'
	
	taxo <- data.frame(kingdom=sapply(tax.tmp.ls, function(x){x[1]}),
					   phylum=sapply(tax.tmp.ls, function(x){x[2]}),
					   class=sapply(tax.tmp.ls, function(x){x[3]}),
					   order=sapply(tax.tmp.ls, function(x){x[4]}),
					   family=sapply(tax.tmp.ls, function(x){x[5]}),
					   genus=sapply(tax.tmp.ls, function(x){x[6]}))
	
	taxo$kingdom <- as.character(taxo$kingdom)
	taxo$phylum <- as.character(taxo$phylum)
	taxo$class <- as.character(taxo$class)
	taxo$order <- as.character(taxo$order)
	taxo$family <- as.character(taxo$family)
	taxo$genus <- as.character(taxo$genus)
	
	for (i in 1:ncol(taxo)){
		taxo[which(is.na(taxo[, i])), i] <- '' 
		}
	
	# taxo.all <- taxo # save big one
	taxo$abundance <- colSums(otu.table)
	row.names(taxo) <- colnames(otu.table)
	
	invisible(taxo)
}


###########  function to bring in QIIME classic OTU table. 
QiimeIn <- function ( file='') {
  tmp.table <- read.table ( file, sep='\t', head=TRUE, row.names=1, comment.char='%', skip=1 )
  tax.numb <- dim(tmp.table)[2]  ## last col is taxon names
  taxa.names <- as.character(tmp.table[,tax.numb])  ## save taxon names
  tmp.to.t <- tmp.table[,1:tax.numb-1]  ## trim names from table
  name <- t(tmp.to.t)  ## transpose after trimming names
  
  preview <- name[1:5,1:5]
  samples <- row.names(name)
  taxa.tmp <- cbind(colnames(name),taxa.names)
  taxa <- data.frame(taxa.tmp)
  names(taxa)[1] <- 'qiime.id'
  row.names(taxa) <- taxa.tmp[, 1]
  dimensions <- dim(name)

  table.info <- list( preview, dimensions, samples, taxa, name )
  names(table.info)[[1]] <- 'Preview'
  names(table.info)[[2]] <- 'Dimensions'
  names(table.info)[[3]] <- 'Samples'
  names(table.info)[[4]] <- 'Taxa'
  names(table.info)[[5]] <- 'Table'
  
  invisible(table.info)
  
}