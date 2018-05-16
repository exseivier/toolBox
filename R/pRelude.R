library(S4Vectors)

#	OBJECTS

#	METHODS

setGeneric("discriminateTargets", function(tar) standarGeneric("discriminateTargets"))
setMethod("discriminateTargets", signature("SimpleList"),
	function(tar) {
		# /* 
		#    It will generates the list of genes that are in common in 2
		#    datasets. Dataset comes from targets object created by selectTargets
		#    function
		# */
		# /* tar
		#    Object created by selectTargets function implemented in pRelude.R
		# */
		tags <- colnames(tar[[1]][["sizes_table"]])[-length(colnames(tar[[1]][["sizes_table"]]))]
		# /*
		#    How many combinations of length 2 do you can do with length(tags) elements?
		# */
		combinations <- comb(2, length(tags))
		results <- SimpleList()
		# // Iterating over the total number of posible combinations of length 2
		for(i in 1:length(combinations[,1])){
			paths <- names(tar)
		# // Iterating over paths stored in tar
			for(p in paths) {
				files <- names(tar[[p]])[-length(names(tar[[p]]))]
				genes_basket <- c()
				out_tag <- paste(tags[combinations[i,1]], tags[combinations[i,2]], sep="")
		# // Iterating over files for ever path sotred in tar
				for(f in files) {
					tag <- tags[combinations[i,1]]
					genes <- tar[[p]][[f]][[tag]]
					tag <- tags[combinations[i,2]]
					genes <- c(genes, tar[[p]][[f]][[tag]])
					genes <- eliminateHeavyDuplicates(names=genes, times=2)
					genes <- gsub("\\.[[:alpha:]]$", "", genes)
					genes <- genes[duplicatedWithMoreThan(names=genes, times=1, drop=FALSE)]
					genes <- unique(sort(genes))
					genes_basket <- c(genes_basket, genes)
				}
				results[[out_tag]] <- t(table(genes_basket))
				results[[paste(out_tag, "_CRgenes", sep="")]] <- unique(sort(genes_basket))
				# /* TODO LIST
				#    Select the Uncoupled regulated genes from L, S and X
				#    in combination of length 2.
				#    Select the coupled and uncoupled regulated genes for
				#    the LSX set
				# */
			}
		}
		results
})

setGeneric("getLists", function(fasta) standardGeneric("getLists"))
setMethod("getLists", signature("character"),
	function(fasta) {
		dna <- readDNAStringSet(fasta, format="fasta")
		genes <- names(dna)
		# /* 
		#    Fasta sequences were removed
		#    to make memory storage efficient.
		# */
		rm(dna)
		# /* extractTags comes from helpe.R
		tags <- extractTags(genes)
		# /*
		#    Calculating homologues and orthlogues genes
		# */
		genes <- unlist(lapply(genes, function(x) strsplit(x, split="\\|")[[1]][1]))
		genes <- genes[grepl("\\.[[:alpha:]]$", genes)]
		alltags <- gsub(".*\\.", "", genes)
		homologs <- SimpleList()
		homologs@listData <- as.list(table(alltags)) # */* Select total L, S and X genes to select singletons genes.
		homologs <- combinations(tags, genes, homologs)
		# /*
		#    select total genes total genes and singletons
		# */
		homologs
})

setGeneric("selectTargets", function(paths, percentage, notags) standardGeneric("selectTargets"))
setMethod("selectTargets", signature("character", "numeric", "numeric"),
	function(paths, percentage, notags) {
		# /* paths.
		#    It is a character vector with the complete paths were *.data files are.
		# */
		# /* percentage.
		#    It is a numeric vector of length 1 which tells the percentage of
		#    the top genes will be selected.
		# */
		# /* notags
		#    It is a numeric vector of length 1 which tells the number of tags included in tar object
		# */
		# /* What it does.
		#    Takes *.data files for each microRNA family.
		#    Select genes by tags.
		#    Sort them by context score from negative to positive.
		#    Takes the $(pecrentage)% of top genes for each subset of genes.
		#    Calculates the vector length of each set of genes, and the total genes
		#    for each tag.
		#    Finally, all results are placed in a SimpleList object.
		# */
		result <- SimpleList()
		for(path in paths) {
			result[[basename(path)]] <- SimpleList()
			lof <- list.files(path=path, pattern="\\.data$")
			min <- 0
			max <- length(lof)
			pb <- txtProgressBar(min=min, max=max, style=3)
			for(f in lof) {
				result[[basename(path)]][[f]] <- SimpleList()
				table <- read.delim(paste(path, f, sep=""), header=FALSE, sep="\t", as.is=TRUE)
				gene_tags <- table[,1]
				gene_tags <- gene_tags[grepl("\\.[[:alnum:]]$", gene_tags)]
				gene_tags <- gsub(".*\\.", "", gene_tags)
				gene_tags <- names(sort(table(gene_tags), decreasing=TRUE))
				for(tag in gene_tags){
					sub_table <- table[grepl(paste("\\.", tag, "$", sep=""), table[,1]),]
					sub_table <- sub_table[order(sub_table[,2], decreasing=FALSE),]
					len_subtable <- length(sub_table[,1])
					sub_table <- sub_table[,1]
					sub_table <- sub_table[1:round(len_subtable*(percentage/100))]
					result[[basename(path)]][[f]][[tag]] <- sub_table
				}
				min <- min + 1
				setTxtProgressBar(pb, min)
			}
		}
		result <- intersect.targets(tar=result, notags=notags)
		result
})

##################################
##################################
##################################


setGeneric("getLists.deprecated", function(total_genes, genes, stage, homoeologos_fasta) standardGeneric("getLists.deprecated"))
setMethod("getLists.deprecated", signature("character", "character", "character", "character"),
	function(total_genes, genes, stage, homoeologos_fasta) {
		list <- SimpleList()
		if (stage == "fu"){
			lof <- list.files(path=genes, pattern="\\.data\\.fu$")
		}
		else if (stage == "fd") {
			lof <- list.files(path=genes, pattern="\\.data\\.fd$")
		}
		else {
			stop(paste(stage, " stage is not understood!", sep=""))
		}
		CRH <- c()
		URH <- c()
		NRH <- c()
		SIN <- c()
		gene_tags <- gsub(".*\\.", "", total_genes)
		gene_tags <- rownames(sort(table(gene_tags), decreasing=TRUE))
		gene_tags <- gene_tags[1:2]
		gene_tags <- sort(gene_tags)
		tgenes <- total_genes[grepl(paste("\\.[", gene_tags[1], gene_tags[2], "]$", sep=""), total_genes)]
		tgenes <- gsub(paste("\\.[", gene_tags[1], gene_tags[2], "]$", sep=""), "", tgenes)
		tgenes <- sort(tgenes)
		tgenes <- eliminateHeavyDuplicates(tgenes, 2)
		homoeologos <- tgenes[duplicated(tgenes)]
		singletons <- tgenes[!(tgenes %in% homoeologos)]
		singletons <- unique(sort(singletons))
		for (f in lof){
			tmp_fu <- read.delim(paste(genes, f, sep=""), header=FALSE, sep="\t", as.is=TRUE)
			tmp_fu <- tmp_fu[,1]
			tmp_fu <- tmp_fu[grepl(paste("\\.[", gene_tags[1], gene_tags[2], "]$", sep=""), tmp_fu)]
			tmp_fu <- gsub(paste("\\.[", gene_tags[1], gene_tags[2], "]$", sep=""), "", tmp_fu)
			tmp_fu <- sort(tmp_fu)
			tmp_fu <- eliminateHeavyDuplicates(tmp_fu, 2)
			tmp_CRH <- tmp_fu[duplicated(tmp_fu)]
			tmp_URH <- tmp_fu[!(tmp_fu %in% tmp_CRH) & tmp_fu %in% homoeologos]
			tmp_SIN <- tmp_fu[!(tmp_fu %in% tmp_CRH) & tmp_fu %in% singletons]
			CRH <- c(CRH, tmp_CRH)
			URH <- c(URH, tmp_URH)
			SIN <- c(SIN, tmp_SIN)
		}
		URH <- URH[!(URH %in% CRH)]
		tCRH <- table(CRH)
		tURH <- table(URH)
		tSIN <- table(SIN)
		CRH <- unique(sort(CRH))
		URH <- unique(sort(URH))
		NRH <- homoeologos[!(homoeologos %in% CRH) & !(homoeologos %in% URH)]
		# Preparing homoeolog genes with reported UTR sequence to filtering NHR. We are extracting the names from homoeologs fasta sequence
		# Be carefully to place the gene name at the begining of the header.
		homoeogenesWUTRs <- readDNAStringSet(homoeologos_fasta, format="fasta")
		homoeogenesWUTRs <- names(homoeogenesWUTRs)
		homoeogenesWUTRs <- lapply(homoeogenesWUTRs, function(x) strsplit(x=x, split="\\|")[[1]])
		items <- c()
		for(item in homoeogenesWUTRs){
			items <- c(items, item[1])
		}
		homoeogenesWUTRs <- items
		homoeogenesWUTRs <- gsub(paste("\\.[", gene_tags[1], gene_tags[2], "]$", sep=""), "", homoeogenesWUTRs)
		homoeogenesWUTRs <- eliminateHeavyDuplicates(homoeogenesWUTRs, 2)
		homoeogenesWUTRs <- homoeogenesWUTRs[duplicated(homoeogenesWUTRs)]
		NRH <- NRH[NRH %in% homoeogenesWUTRs]
		NRH <- eliminateHeavyDuplicates(NRH, 1)
		URH <- URH[URH %in% homoeogenesWUTRs]
		URH <- eliminateHeavyDuplicates(URH, 1)
		list[["total_genes"]] <- total_genes
		list[["homoeologos"]] <- homoeologos
		list[["singletons"]] <- singletons
		list[["CRH"]] <- CRH
		list[["URH"]] <- URH
		list[["NRH"]] <- NRH
		list[["tCRH"]] <- tCRH
		list[["tURH"]] <- tURH
		list[["tSIN"]] <- tSIN
		print("debbuging2")
		list
})

setGeneric("loadTargets", function(path, lof, collapse=TRUE) standarGeneric("loadTargets"))
setMethod("loadTargets", signature("character", "character"),
	# /*
	#    Load the names of the target genes from each *.data.fu or *.data.fd file.
	#    Those files are the output of the target scan analysis and *.fu are the
	#    best prediction targets and *.fd are the worst prediction targets.
	# */
	# /* path
	#    Character vector with the name of the path where the *.data.f* are.
	# */
	# /* lof
	#    Character vector with the names of the *.data.f* files
	# */
	# /* collapse
	#    Logical vector. If TRUE gene names will be filtered and untagged, sorted
	#    and unique-ed
	function(path, lof, collapse=TRUE) {
		if(length(path) == 0 | length(lof) == 0){
			stop("Path or list of files were not specified")
		}
		results <- SimpleList()
		for (f in lof) {
			tmp_data <- read.delim(paste(path, f, sep=""), header=FALSE, sep="\t", as.is=TRUE)[,1]
			tmp_data <- gsub("\\|.*$", "", tmp_data)
			if (collapse) {
				tmp_data <- gsub("\\.[LSPX]$", "", tmp_data)
				tmp_data <- tmp_data[tmp_data %!in% "unnamed"]
				tmp_data <- unique(sort(tmp_data))
			}
			tmp_data <- tmp_data[tmp_data %!in% "unnamed"]
					results[[ strsplit(f, split="\\.")[[1]][1]  ]] <- tmp_data
		}
		results
	
})

setGeneric("filter.bad.targs", function(targs, badtargs) standarGeneric("filter.bad.targs"))
setMethod("filter.bad.targs", signature("SimpleList", "SimpleList"),
	# /*
	#
	# */
	function(targs, badtargs) {
		realbadtargs <- SimpleList()
		alltargs <- c()
		for (name in names(targs)) {
			alltargs <- c(alltargs, targs[[name]])
		}
		for (name in names(badtargs)) {
			realbadtargs[[name]] <- badtargs[[name]][badtargs[[name]] %!in% alltargs]
		}
		realbadtargs
})
