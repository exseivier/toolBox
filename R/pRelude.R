library(S4Vectors)

#	OBJECTS

#	METHODS

setGeneric("getLists", function(total_genes, genes, stage) standardGeneric("getLists"))
setMethod("getLists", signature("character", "character", "character"),
	function(total_genes, genes, stage) {
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
		tgenes <- gsub("\\.[LS]$", "", total_genes)
		tgenes <- sort(tgenes)
		homoeologos <- unique(sort(tgenes[duplicated(tgenes)]))
		singletons <- tgenes[!(tgenes %in% homoeologos)]
		singletons <- unique(sort(singletons))
		for (f in lof){
			tmp_fu <- read.delim(paste(genes, f, sep=""), header=FALSE, sep="\t", as.is=TRUE)
			tmp_fu <- tmp_fu[,1]
			tmp_fu <- gsub("\\.[LS]$", "", tmp_fu)
			tmp_fu <- sort(tmp_fu)
			tmp_CRH <- unique(sort(tmp_fu[duplicated(tmp_fu)]))
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
		list[["total_genes"]] <- total_genes
		list[["homoeologos"]] <- homoeologos
		list[["singletons"]] <- singletons
		list[["CRH"]] <- CRH
		list[["URH"]] <- URH
		list[["NRH"]] <- NRH
		list[["tCRH"]] <- tCRH
		list[["tURH"]] <- tURH
		list[["tSIN"]] <- tSIN
		list
})
