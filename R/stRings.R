#	This library is intended to measure properties of strings and to manipulate them

#	OBJECTS


#	METHODS
setGeneric("getPairLenDiff", function(path, file1, file2) standardGeneric("getPairLenDiff"))
setMethod("getPairLenDiff", signature("character", "character", "character"),
	function(path, file1, file2) {
		seqs <- SimpleList()
		seq1 <- readDNAStringSet(paste(path, file1, sep="/"), format="fasta")
		seq1 <- seq1[order(names(seq1))]
		seq2 <- readDNAStringSet(paste(path, file2, sep="/"), format="fasta")
		seq2 <- seq2[order(names(seq2))]
		seqs[[file1]] <- seq1
		seqs[[file2]] <- seq2
		if (length(seqs@listData[[file1]]) != length(seqs@listData[[file2]])) {
			stop("Number of sequences of file1 does not match the number of sequences of file2")
		}
		cat("\nPreparing orthogroup genes\n")
		min <- 0
		max <- length(seqs@listData[[1]])
		pb <- txtProgressBar(min=min, max=max, style=3)

		gNames <- names(seqs@listData[[1]])
		gNames <- sapply(lapply(strsplit(gNames, ""), rev), paste, collapse="")
		gNames <- gsub("^.*\\|[[:alnum:]]*\\.", "", gNames)
		gNames <- sapply(lapply(strsplit(gNames, ""), rev), paste, collapse="")
		orthoGroups <- SimpleList()
		for (i in 1:length(seqs@listData[[1]])) {
			seqList <- c()
			names <- c()
			for (j in 1:length(seqs@listData)) {
				seqList <- c(seqList, toString(seqs@listData[[j]][i]))
				names <- c(names, names(seqs@listData[[j]][i]))
			}
			names(seqList) <- names
			orthoGroups@listData[[gNames[i]]] <- DNAStringSet(seqList, use.names=TRUE)
			setTxtProgressBar(pb, i)
		}

		cat("\nCalculating sequence-length difference between orthologs\nFile1 Vs File2\n")
		min <- 0
		max <- length(orthoGroups)
		pb <- txtProgressBar(min=min, max=max, style=3)
		gNames <- names(orthoGroups)
		lenDiff <- SimpleList()
		for (i in 1:length(orthoGroups)) {
			for (j in 1:(length(orthoGroups[[i]])-1)) {
				lenDiff@listData[[gNames[i]]] <- width(orthoGroups@listData[[i]][j]) - width(orthoGroups@listData[[i]][j+1])
			}
			setTxtProgressBar(pb, i)
		}
		lenDiff
})

setMethod("getPairLenDiff", signature("SimpleList", "character", "character"),
	function(path, file1, file2) {
		
})

