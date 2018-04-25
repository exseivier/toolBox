# OBJECTS

# METHODS

setGeneric("create.intersectIO", function(paths, tag, homoeologos) standardGeneric("create.intersectIO"))
setMethod("create.intersectIO", signature("character", "character", "character"),
	function(paths, tag, homoeologos) {
		# homoeologos is a character vector with the names of homoeologs genes
		# tag is a character vector with the tag of .data out files of target scan
		# paths is a character vector with the paths where the *.data files are
		results <- SimpleList()
		for (i in 1:length(paths)) {
			lof <- list.files(path=paths[i], pattern=paste("\\.data\\.", tag, "$", sep="")) # lof: list of files
			path <- basename(paths[i])
			results[[path]] <- SimpleList()
			for (f in lof) { # For each file in list of files
				results[[path]][[f]] <- SimpleList()
				genes <- read.delim(paste(paths[i], f, sep=""), header=FALSE, sep="\t", as.is=TRUE)
				genes <- genes[,1]
				gene_tags <- gsub(".*\\.", "", genes)
				gene_tags <- rownames(sort(table(gene_tags), decreasing=TRUE))
				gene_tags <- gene_tags[1:2]
				gene_tags <- sort(gene_tags)
				genes <- genes[gsub(paste("\\.[", paste(gene_tags[1], gene_tags[2], sep=""), "]$", sep=""), "", genes) %in% homoeologos]
				results[[path]][[f]][[gene_tags[1]]] <- genes[grepl(paste("\\.", gene_tags[1], "$", sep=""), genes)]
				results[[path]][[f]][[gene_tags[2]]] <- genes[grepl(paste("\\.", gene_tags[2], "$", sep=""), genes)]
				genes <- gsub(paste("\\.[", paste(gene_tags[1], gene_tags[2], sep=""), "]$", sep=""), "", genes)
				genes <- sort(genes)
				genes <- genes[duplicated(genes)]
				results[[path]][[f]][[paste(gene_tags[1], gene_tags[2], sep="")]] <- genes
			}
		}
		mirnas_families <- names(results[[1]])

		for (family in mirnas_families) {
			genes <- c()
			for (i in 1:length(results)) {
				genes <- c(genes, results[[i]][[family]][[3]])
			}
			genes <- sort(genes)
#			genes <- genes[duplicated(genes)]
#			genes <- sort(genes)
#			genes <- genes[duplicated(genes)]
#			genes <- unique(genes)
			genes <- eliminateHeavyDuplicates(genes, length(results))
			genes <- genes[duplicatedWithMoreThan(genes, length(results)-1)]
			results[["full_intersection"]][[family]] <- genes
		}
		results
})

setGeneric("genrich", function(gList, interIO, tag) standardGeneric("genrich"))
setMethod("genrich", signature("SimpleList", "SimpleList", "character"),
	function(gList, interIO, tag) {
		#	gList is an object created by getLists function.
		#	inerIO is an object created by create.intersectIO function.
		#	tag is a character vector with the name of the dataset you want to use which comes from interIO object.
		genetags <- gList[["total_genes"]]
		genetags <- gsub(".*\\.", "", genetags)
		genetags <- sort(table(genetags), decreasing=TRUE)
		genetags <- names(genetags)[1:2]
		genetags <- sort(genetags)
		balls_in_urn <- c(gList[["CRH"]], gList[["URH"]], gList[["NRH"]])
		balls_in_urn <- 2*length(balls_in_urn)
		interIO <- interIO[[tag]]
		m <- balls_in_urn / 2	# White balls in urn
		n <- m					# Black balls in urn
		pvalues <- SimpleList()
		enriched_genes <- c()
		for(f in names(interIO)){
			pvalues[[f]] <- SimpleList()
			k <- length(interIO[[f]][[genetags[1]]]) + length(interIO[[f]][[genetags[2]]])	# Balls in sample
			x <- length(interIO[[f]][[genetags[1]]])
			pvalue <- -phyper(m=m, n=n, k=k, q=x, log.p=TRUE, lower.tail=FALSE)
			pvalues[[f]][[genetags[1]]] <- pvalue
			x <- length(interIO[[f]][[genetags[2]]])
			pvalue <- -phyper(m=m, n=n, k=k, q=x, log.p=TRUE, lower.tail=FALSE)
			pvalues[[f]][[genetags[2]]] <- pvalue
			if (pvalues[[f]][[genetags[1]]] > 2){
				pvalues[[f]][["enriched_genes"]] <- genetags[1]
			}
			else if (pvalues[[f]][[genetags[2]]] > 2) {
				pvalues[[f]][["enriched_genes"]] <- genetags[2]
			}
			else {
				pvalues[[f]][["enriched_genes"]] <- "NONE"
			}
		}
		data <- c()
		for(i in 1:length(pvalues)) {
			row_data <- c()
			for(j in 1:length(pvalues[[i]])){
				row_data <- c(row_data, pvalues[[i]][[j]])
			}
			data <- rbind(data, row_data)
		}
		rownames(data) <- names(pvalues)
		data
})
