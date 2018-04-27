
# FUNCTIONS
eliminateHeavyDuplicates <- function(names, times) {
	tmp_names <- sort(names)
	substracting_names <- sort(names)
	for(i in 1:times) {
		substracting_names <- substracting_names[duplicated(substracting_names)]
	}
	tmp_names <- tmp_names[!tmp_names %in% substracting_names]
	tmp_names
}

duplicatedWithMoreThan <- function(names, times, drop=FALSE) {
	tmp_names <- sort(names)
	for(i in 1:times) {
		tmp_names <- tmp_names[duplicated(tmp_names)]
	}
	if (drop == FALSE) {
		names %in% tmp_names
	}
	else {
		tmp_names
	}
}

testSecondDuplicate<- function(vec){
	saw <- c()
	bool <- c()
	for(i in 1:length(vec)) {
		if (vec[i] %in% saw){
			bool <- c(bool, TRUE)
		}
		else {
			bool <- c(bool, FALSE)
			saw <- c(saw, vec[i])
		}
	}
	bool
}

comb <- function(size, universe) {
	numseq <- as.numeric(paste(rep(1, size), collapse="")):as.numeric(paste(rep(universe, size), collapse=""))
	numseq <- as.character(numseq)
	m <- matrix(unlist(lapply(numseq, function(x) as.numeric(strsplit(x, split="")[[1]]))), ncol=size, byrow=TRUE)
	m <- m[apply(m, 1, function(x) !any(x > universe | x == 0)),]
	m <- m[apply(apply(m, 1, function(x) duplicatedWithMoreThan(x, times=1, drop=FALSE)), 2, sum) == 0,]
	for(i in 1:(size-1)) {
		print(i)
		m <- m[m[,i+1] > m[,i],]
	}
	m
}

intersect.targets <- function(tar, notags=3) {
	# /* tar is an object created by selectTargets function from pRelude.R */
	# /* notags is a numeric vector of length 1 which tells the number of tags included in tar object */
	if(class(tar) != "SimpleList" & class(notags) != "numeric") {
		stop("The argumet is a bad data structure")
	}
	sizes <- SimpleList()

	for (i in 1:length(tar)){						#	// ITERATING ALONG PATHS
		names <- names(tar[[i]])
		path <- names(tar)[i]
		for(f in 1:length(tar[[path]])){			#	// ITERATING ALONG FILENAMES
			file <- names(tar[[path]])[f]
			for(t in 1:length(tar[[path]][[file]])){		#	// ITERATING ALONG TAGS
				tag <- names(tar[[path]][[file]])[t]
				if (tag %in% names(sizes)){
					sizes[[ tag ]] <- c(sizes[[ tag ]], length(tar[[path]][[file]][[tag]]))
				}
				else {
					sizes[[ tag ]] <- c(length(tar[[path]][[file]][[tag]]))
				}
			}
		}
		# Making sure that elements of the list have the same length
		sizes <- sizes[names(sizes) %in% names(sapply(sizes,length)[duplicatedWithMoreThan(names=as.character(sapply(sizes, length)), times=notags-1, drop=FALSE)])]
		names_sizes <- names(sizes)
		sizes <- matrix(unlist(sizes@listData), ncol=length(sizes), byrow=FALSE)
		colnames(sizes) <- names_sizes
		rownames(sizes) <- names
		sizes <- cbind(sizes, total_sum=apply(sizes, 1, sum))
		tar[[path]][["sizes_table"]] <- SimpleList()
		tar[[path]][["sizes_table"]] <- sizes

		# /* TODO list
		#    For every element in tar object perform the calculus of
		#    intersection percentage in all posible combinations of length
		#    2 to n-1, where n is the length of the tags of every gene set.
		#    Add the results of this calculus to tar object
		# */
	
	}
	tar
}

extractTags <- function(names) {
	names <- unlist(lapply(names, function(x) strsplit(x, split="\\|")[[1]][1]))
	tags <- names[grepl("\\.[[:alpha:]]$", names)]
	tags <- gsub(".*\\.", "", tags)
	tags <- table(tags)
	tags <- names(tags)
	tags
}

combinations.internals <- function(tmp_tags, tmp_genes) {
	tmp_genes <- tmp_genes[ grepl( paste( "\\.[", paste(tmp_tags, collapse=""), "]$", sep="" ), tmp_genes ) ]
	tmp_genes <- gsub(paste("\\.[", paste(tmp_tags, collapse=""), "]$", sep=""), "", tmp_genes)
	tmp_genes <- tmp_genes[duplicatedWithMoreThan(tmp_genes, times=(length(tmp_tags)-1), drop=FALSE)]
	tmp_genes <- eliminateHeavyDuplicates(tmp_genes, times=length(tmp_tags))
	tmp_genes <- unique(sort(tmp_genes))
	tmp_genes
}

combinations <- function(tags, genes, homologs) {
	tmp_tags <- tags
	for (i in 1:(length(tags)-1)) {
		tmp_genes <- genes
		homologs[[ paste(tmp_tags, collapse="") ]] <- combinations.internals(tmp_tags, tmp_genes)
		tmp_tags <- tmp_tags[1:(length(tmp_tags)-1)]
	}
	com <- comb(2, length(tags))
	for(r in 2:length(com[,1])){
		row <- com[r,]
		tmp_tags <- c(tags[com[r,1]], tags[com[r,2]])
		tmp_genes <- genes
		homologs[[ paste(tmp_tags, collapse="") ]] <- combinations.internals(tmp_tags, tmp_genes)
	}
	homologs
}
