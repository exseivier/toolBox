library(S4Vectors)

#	OBJECTS

#	METHODS

setGeneric("getLists", function(total_genes, genes, stage, homoeologos_fasta) standardGeneric("getLists"))
setMethod("getLists", signature("character", "character", "character", "character"),
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
		tgenes <- total_genes[grepl("\\.[LS]$", total_genes)]
		tgenes <- gsub("\\.[LS]$", "", tgenes)
		tgenes <- sort(tgenes)
		tgenes <- eliminateHeavyDuplicates(tgenes, 2)
		homoeologos <- tgenes[duplicated(tgenes)]
		singletons <- tgenes[!(tgenes %in% homoeologos)]
		singletons <- unique(sort(singletons))
		for (f in lof){
			tmp_fu <- read.delim(paste(genes, f, sep=""), header=FALSE, sep="\t", as.is=TRUE)
			tmp_fu <- tmp_fu[,1]
			tmp_fu <- tmp_fu[grepl("\\.[LS]$", tmp_fu)]
			tmp_fu <- gsub("\\.[LS]$", "", tmp_fu)
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
		homoeogenesWUTRs <- gsub("\\.[LS]$", "", homoeogenesWUTRs)
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
		print("debbuging1")
		list
})


