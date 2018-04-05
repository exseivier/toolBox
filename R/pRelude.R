#	OBJECTS

#	METHODS

setGeneric("getLists", function(total_genes, genes) standardGeneric("getLists"))
setMethod("getLists", signature("character", "character"),
	function(total_genes, genes) {
		lists <- SimpleList()
		total_genes <- gsub("\\.[LS]$", "", total_genes)
		genes <- gsub("\\.[LS]$", "", genes)
		homoeologos <- total_genes[unique(duplicated(sort(total_genes)))]
		singletons <- total_genes[!(total_genes %in% homoeologos)]
		CRH <- genes[unique(duplicated(sort(genes)))]
		CURH <- genes[!(genes %in% singletons)]
		URH <- CURH[!(CURH %in% CRH)]
		lists[["homoeologos"]] <- homoeologos
		lists[["singletons"]] <- singletons
		lists[["CRH"]] <- CRH
		lists[["CURH"]] <- CURH
		lists[["URH"]] <- URH
		lists[["total_genes"]] <- total_genes
		lists[["genes"]] <- genes
		lists[["tbl_genes"]] <- as.matrix(table(genes))
		colnames(lists[["tbl_genes"]]) <- "No. microRNAs families"
		lists
})


