# OBJECTS

setClass("intersectio", representation(listData="list", names="character"), prototype(listData=list(), names=NA_character_),
	validity=function(object){
		if ( length(object@listData) == 0) {
			stop("listData is empty")
		}
		if (length(object@names) <= 1) {
			if(is.na(object@names)) {
				stop("Gene names vector is empty")
			}
		}
#		for (item in object@listData) {
#			if (length(object@names) != length(item[[1]][,1])) {
#				print(str(item))
#				print(paste("Gene-names vector length: ", length(object@names), sep=""))
#				stop("Gene-names vector length does not match with the column length in data frame!")
#			}
#		}
})


# METHODS

# Matching methods

setGeneric("match.rho2time", function(mat1, mat2, scaled_time_pos_table, colname, cmd) standardGeneric("match.rho2time"))
setMethod("match.rho2time", signature(mat1="matrix", mat2="matrix", scaled_time_pos_table="data.frame", colname="character", cmd="character"),
	function(mat1, mat2, scaled_time_pos_table, colname, cmd) {
		result <- c()
		if (is.null(rownames(mat1)) || is.null(rownames(mat2))) {
			stop("No rownames in one or both matrices mat1 or mat2")
		}
		rhoQ <- mat1[,colname] / mat2[,colname]
		rhoD <- mat1[,colname] - mat2[,colname]
		rho <- cbind(rhoQ, rhoD, mat1[,colname], mat2[,colname])
		time_pos <- c()
		for (i in 1:length(rho[,1])) {
			time_pos <- c(time_pos, scaled_time_pos_table[scaled_time_pos_table[,2] == rownames(rho)[i], 4])
		}
		rho <- cbind(rho, scaled_time_pos_table)
		if (cmd == "plot") {
			boxplot(rho[,1] ~ rho[,5], ylab="rho Quotient", xlab="Scaled time position")
			atripchart(rho[,1] ~ rho[,5], col="red", pch=20, method="jitter", add=T, vertical=T)
		}
		else {
			rho
		}
})

setGeneric("match.genes2UTRlen", function(intergenes, genes2UTRlen) standardGeneric("match.genes2UTRlen"))
setMethod("match.genes2UTRlen", signature("matrix", "data.frame"),
	function(intergenes, genes2UTRlen) {
	genes2UTRlen <- genes2UTRlen[gsub("\\.[LS]$", "", genes2UTRlen[,1]) %in% rownames(intergenes),]
	result <- c()
	for(i in 1:length(genes2UTRlen[,1])) {
		gene_name <- as.character(genes2UTRlen[i,1])
		seq_len <- as.numeric(genes2UTRlen[i,2])
		freq <- as.numeric(intergenes[rownames(intergenes) == gsub("\\.[LS]$", "", gene_name),1])
		result <- rbind(result, c(gene_name, seq_len, freq))
	}
	row_names <- result[,1]
	result <- result[,2:3]
	result <- apply(result, 2, as.character)
	result <- apply(result, 2, as.numeric)
	rownames(result) <- row_names
	colnames(result) <- c("seq_len", "freq")
	result
})


# Methods for creating objects

setGeneric("create.intersectio.Object.2", function(paths) standardGeneric("create.intersectio.Object.2"))
setMethod("create.intersectio.Object.2", signature("character"),
	function(paths){
		DATA <- list()
		for(path in paths){
			files <- list.files(path=path, pattern="*.data$")
			COUNTS_BEST <- list()
			COUNTS_WORST <- list()
			TARGETS <- list()
			INTERSECT <- list()
			print(paste("Processing ", path, "...", sep=""))
			min <- 0
			pb <- txtProgressBar(min=min, max=length(files), style=3)
			for (file in files) {
				all_genes <- read.delim(paste(path, file, sep="/"), header=F, ) # Reading all genes and context scores
				percent_25 <- round(0.25 * length(all_genes[,1]), 0)
				best_pred <- all_genes[(1:percent_25),] # Taking the best predicted targets
				worst_pred <- all_genes[(length(all_genes[,1]) - percent_25):length(all_genes[,1]),] # Taking the worst predicted targets
				command <- paste(names(table(gsub(".+\\.", "", all_genes[,1]))), collapse="")
				best_pred <- split(best_pred, gsub(".+\\.", "", best_pred[,1])) # splitting Best predicted genes L, S or X. It depends on the tag attached at the end of the gene name
				worst_pred <- split(worst_pred, gsub(".+\\.", "", worst_pred[,1])) # Splitting Worst predicted genes in the same way like the Best predicted genes
				

				for (i in 1:length(best_pred)) {	# Counting best targets
					best_pred[[i]][,1] <- gsub(paste("\\.[", command, "]$", sep=""), "", best_pred[[i]][,1])
					TARGETS[[names(best_pred)[i]]] <- best_pred[[i]][,1]
				}
				COUNTS_BEST[[file]] <- list(TARGETS=TARGETS)
				for (i in 1:length(worst_pred)) {	# Counting worst targets
					worst_pred[[i]][,1] <- gsub(paste("\\.[", command, "]$", sep=""), "", worst_pred[[i]][,1])
					TARGETS[[names(worst_pred)[i]]] <- worst_pred[[i]][,1]
				}
				COUNTS_WORST[[file]] <- list(TARGETS=TARGETS)

				for (i in 1:(length(best_pred)-1)) {	# Counting Intersections for best predictions
					tag <- paste(names(best_pred)[i],  names(best_pred)[i+1], collapse="")
					INTERSECT[[tag]] <- intersect(best_pred[[i]][,1], best_pred[[i+1]][,1])
				}
				COUNTS_BEST[[file]] <- c(COUNTS_BEST[[file]], list(INTERSECT=INTERSECT))
				
				for (i in 1:(length(worst_pred)-1)) {
					tag <- paste(names(worst_pred)[i], names(worst_pred)[i+1], collapse="")
					INTERSECT[[tag]] <- intersect(worst_pred[[i]][,1], worst_pred[[i+1]][,1])
				}
				COUNTS_WORST[[file]] <- c(COUNTS_WORST[[file]], list(INTERSECT=INTERSECT))
				min <- min + 1
				setTxtProgressBar(pb, min)

			}
			DATA[[paste("BEST_", command, sep="")]] <- COUNTS_BEST
			DATA[[paste("WORST_", command, sep="")]] <- COUNTS_WORST
		}
		result <- new("intersectio", listData=DATA, names=files)
		result
})


setGeneric("create.intersectio.Object", function(paths, commands, tag) standardGeneric("create.intersectio.Object"))
setMethod("create.intersectio.Object", signature("character", "character", "character"),
	function(paths, commands, tag){
		if(length(paths) != length(commands)) {
			stop("paths length is different to commands length")
		}
		i <- 1
		data <- list()
		for(path in paths) {
			files.fu <- list.files(path=path, pattern=paste("^[ACGTU].*\\.", tag, "$", sep=""))
			#files.fd <- list.files(path=path, pattern="^[ACGTU].*\\.fd$")
			command <- commands[i]
			i <- i + 1
			x <- c()
			y <- c()
			z <- c() # for the intersection
			t <- c() # for total
			intergenes <- list()
			geneList <- list()	### Mod 1 20180122 ###
			counts <- c()
			for (file in files.fu) {
				genes <- read.delim(paste(path, file, sep=""), header=F, row.names=1, sep="\t")
				f_com <- strsplit(command, split="")[[1]][1]
				s_com <- strsplit(command, split="")[[1]][2]
				genes <- rownames(genes)
				geneList[[file]] <- genes	### Mod 1 20180122 ###
				t <- c(t, length(genes))
				x <- c(x, length(genes[grepl(paste("\\.[", f_com, "]$", sep=""), genes)]))
				y <- c(y, length(genes[grepl(paste("\\.[", s_com, "]$", sep=""), genes)]))
				genes <- gsub(paste("\\.[", f_com, s_com, "]$", sep=""), "", genes)
				genes <- sort(genes)
				genes <- genes[duplicated(genes)]
				genes <- unique(genes)
			#	geneList[[file]] <- geneList[[file]][!(geneList[[file]] %in% genes)]	### Mod 1 20180122 ###
				z <- c(z, length(genes))
				intergenes[[file]] <- genes
			}
			counts <- cbind(x, y, z, t)
			rownames(counts) <- files.fu
			colnames(counts) <- c(f_com, s_com, command, "Total")
			data[[command]] <- list(counts=counts, intergenes=intergenes, NoIntergenes=geneList)
		}
		result <- new("intersectio", listData=data, names=files.fu)
		result
})

setGeneric("prepare.intersectio", function(intersectio) standardGeneric("prepare.intersectio"))
setMethod("prepare.intersectio", signature("intersectio"), 
	function(intersectio) {
		if (!("IALL" %in% names(intersectio@listData))) {
			print("This object does not have the IALL slot")
			stop("Prepare it with intersect.intersectio function")
		}
		result <- c()
		rownames <- intersectio@names
		ncol <- (length(intersectio@listData)-1)
		nrow <- length(intersectio@listData[[1]]$counts[,1])
		totals <- matrix(rep(0, ncol * nrow), ncol=ncol, nrow=nrow)
		for (i in 1:(length(intersectio@listData)-1)) {
			table <- c() # Used to store the percentages of intersection and L and S
			tot_pos <- length(colnames(intersectio@listData[[i]]$counts)) # total counts column position
			int_pos <- (length(colnames(intersectio@listData[[i]]$counts)) - 1) # Intersection counts column position
			for (j in 1:length(intersectio@listData[[i]]$counts[,1])) {
				row <- intersectio@listData[[i]]$counts[j,]
				totals[j,i] <- row[tot_pos]
				result_row <- c() # Used to store the percentage results by row
				for (k in 1:(int_pos-1)) {
					result_row <- c(result_row, ((row[k] - row[int_pos]) / row[tot_pos]) * 100) # Calculating percentages for L, S or X
				}
				result_row <- c(result_row, (((row[int_pos] * 2) / row[tot_pos]) * 100)) # Calculating the percentage of intersection
				table <- rbind(table, result_row)
			}
			result <- cbind(result, table)
		}
		# Calculating percentage of LSX intersection
		num_intersect <- as.vector(intersectio@listData$IALL$counts[,1])
		#Checking compatibility between vector and matrix
		if(length(num_intersect) != length(totals[,1])) {
			stop("Length of num_intersect vector does not match with length of totals[,1]")
		}
		percent_num_intersect <- c()
		for (i in 1:length(totals[,1])) {
			sum <- 0
			for (j in 1:length(totals[1,])) {
				sum <- sum + totals[i,j]
			}
			percent_num_intersect <- c(percent_num_intersect, (((length(totals[1,])*num_intersect[i]) / sum) * 100))
		}
		result <- cbind(result, IALL=percent_num_intersect)
		rownames(result) <- rownames
		result
})

setGeneric("intersect.intersectio", function(intersectio) standardGeneric("intersect.intersectio"))
setMethod("intersect.intersectio", signature("intersectio"),
	function(intersectio) {
		data <- list()
		for (item in intersectio@listData) {
			file_names <- names(item$intergenes)
			i <- 1
			for (intergene in item$intergenes) {
				data[[file_names[i]]] <- c(data[[file_names[i]]], intergene)
				i <- i + 1
			}
			names(data) <- file_names
		}
		result <- c()
		intergenes <- list()
		for (id in names(data)) {
			data[[id]] <- sort(data[[id]])
			data[[id]] <- data[[id]][duplicated(data[[id]])]
			data[[id]] <- unique(data[[id]])
			result <- rbind(result, length(data[[id]]))
		}
		rownames(result) <- names(data)
		colnames(result) <- c("IALL")
		intersectio@listData$IALL <- list(counts=result, intergenes=data)
		result <- intersectio
		result
})

# Filtering methods
setGeneric("filter.NoIntergenes.UTRlen", function(intergenes, genes2UTRlen) standardGeneric("filter.NoIntergenes.UTRlen"))
setMethod("filter.NoIntergenes.UTRlen", signature("matrix", "data.frame"),
	function(intergenes, genes2UTRlen){
	genes2UTRlen <- genes2UTRlen[!gsub("\\.[LS]$", "", genes2UTRlen[,1]) %in% rownames(intergenes),]
	genes2UTRlen
})

setGeneric("filter.single_and_NoRAAC.UTRlen", function(intergenes, genes2UTRlen) standardGeneric("filter.single_and_NoRAAC.UTRlen"))
setMethod("filter.single_and_NoRAAC.UTRlen", signature("matrix", "data.frame"),
	function(intergenes, genes2UTRlen) {
	result <- list()
	noIntergenes <- filter.NoIntergenes.UTRlen(intergenes=intergenes, genes2UTRlen=genes2UTRlen)
	noIntergenes <- noIntergenes[,1]
	noIntergenes <- as.vector(noIntergenes)
	noIntergenes <- gsub("\\.[LS]$", "", noIntergenes)
	noIntergenes <- sort(noIntergenes)
	noI_dup <- noIntergenes[duplicated(noIntergenes)]
	noI_sin <- noIntergenes[!duplicated(noIntergenes)]
	genes2UTRlen_dup <- genes2UTRlen[gsub("\\.[LS]$", "", genes2UTRlen[,1]) %in% noI_dup,]
	genes2UTRlen_sin <- genes2UTRlen[gsub("\\.[LS]$", "", genes2UTRlen[,1]) %in% noI_sin,]
	result[["NoRAAC"]] <- genes2UTRlen_dup
	result[["singleton"]] <- genes2UTRlen_sin
	result
})

# Hypothesis testing methods
setGeneric("kruskal.do", function(data_1, data_2) standardGeneric("kruskal.do"))
setMethod("kruskal.do", signature("numeric", "numeric"),
	function(data_1, data_2){
		first <- rep("first", length(data_1))
		second <- rep("second", length(data_2))
		factors <- c(first, second)
		data <- c(data_1, data_2)
		data <- cbind(data, factors)
		data <- as.data.frame(data)
		data$data <- as.numeric(data$data)
		data$factors <- factor(data$factors)
		kruskal <- kruskal.test(data$data, data$factors)
		kruskal
})

setGeneric("GOF.dunif", function(x, y) standardGeneric("GOF.dunif"))
setMethod("GOF.dunif", signature("numeric", "numeric"),
	function(x, y){
		if(length(x) != length(y)) {
			stop("x vector length does not match with y vector length")
		}
		new_y <- y * x
		print(cbind(new_y, y, x))
		result <- cor.test(x, new_y)
		result
})

setGeneric("shift.intergenes", function(sorted_genes, targets, intergenes, width, cmd) standardGeneric("shift.intergenes"))
setMethod("shift.intergenes", signature("character", "character", "character", "numeric", "character"),
	function(sorted_genes, targets, intergenes, width, cmd) {
	tag_interS <- c()
	tag_interL <- c()
	for(i in 1:length(sorted_genes)) {
		tag_interL <- c(tag_interL, ifelse(gsub("\\.[LS]$", "", sorted_genes[i]) %in% intergenes && grepl("\\.L$", sorted_genes[i]), 1, 0))
		tag_interS <- c(tag_interS, ifelse(gsub("\\.[LS]$", "", sorted_genes[i]) %in% intergenes && grepl("\\.S$", sorted_genes[i]), 1, 0))
	}
	total <- length(sorted_genes)
	total_interL <- length(tag_interL[tag_interL == 1]) # M / N White or Black balls in urn
	total_interS <- length(tag_interS[tag_interS == 1]) # M / N White or Black balls in urn
#	total_inter <- total_interL + total_interS # M + N total balls in urn
	logpv.L <- c()
	logpv.S <- c()
	x_axis <- c()
	for (i in 1:(total-width)) {
		sub_interL <- tag_interL[i:(i+width)]
		sub_interL <- length(sub_interL[sub_interL == 1]) # X for L
		sub_interS <- tag_interS[i:(i+width)]
		sub_interS <- length(sub_interS[sub_interS == 1]) # X for S
		sub_inter <- sub_interL + sub_interS # K sampled balls
		enrich_L <- -phyper(sub_interL, m=total_interL, n=total_interS, k=sub_inter, log.p=T, lower.tail=F)
		logpv.L <- c(logpv.L, enrich_L)
		enrich_S <- -phyper(sub_interS, m=total_interS, n=total_interL, k=sub_inter, log.p=T, lower.tail=F)
		logpv.S <- c(logpv.S, enrich_S)
	}
	par(mfrow=c(2,1))
	y_max <- max(max(logpv.L[is.finite(logpv.L)]), max(logpv.S[is.finite(logpv.S)]))
	y_min <- min(min(logpv.L[is.finite(logpv.L)]), min(logpv.S[is.finite(logpv.S)]))
	plot(1:length(sorted_genes), tag_interL, col="blue", type="l", lwd=1)
	lines(1:length(sorted_genes), tag_interS, col="red", lty=1, lwd=1)
	plot((width/2):(total-((width/2)+1)), logpv.L, col="blue", type="l", lwd=1, xlim=c(0, length(sorted_genes)), ylim=c(y_min, y_max))
	lines((width/2):(total-((width/2)+1)), logpv.S, col="red", lty=1, lwd=1)
	abline(h=2, col="black", lwd=2)
})

setGeneric("test.GO", function(gene2cat, all_genes, intersectio, genome_file, method, cmd) standardGeneric("test.GO"))
setMethod("test.GO", signature("data.frame", "character", "intersectio", "character", "character", "character"),
	function(gene2cat, all_genes, intersectio, genome_file, method, cmd) {
		intergenes_for_all_mirnas <- c() # Also known as RAAC genes
		for (item in intersectio@listData$LS$intergenes) {
			intergenes_for_all_mirnas <- c(intergenes_for_all_mirnas, item)
		}
		intergenes_for_all_mirnas <- sort(intergenes_for_all_mirnas)
		intergenes_for_all_mirnas <- unique(intergenes_for_all_mirnas) # RAAC genes: homoeologous regulated as a couple
		print(length(intergenes_for_all_mirnas))

		all_genes <- gsub("\\.[LS]$", "", all_genes)
		all_genes <- sort(all_genes)
		NoRAAC <- all_genes[duplicated(all_genes)]
		NoRAAC <- unique(NoRAAC)
		NoRAAC <- NoRAAC[!(NoRAAC %in% intergenes_for_all_mirnas)] # NoRAAC genes: homoeologous not regulated as a couple

		single <- all_genes[!(all_genes %in% NoRAAC) & !(all_genes %in% intergenes_for_all_mirnas)]
		
		print("Intergenes in single genes?")
		print(sum(intergenes_for_all_mirnas %in% single))
		print("Intergenes in NoRAAC genes?")
		print(sum(intergenes_for_all_mirnas %in% NoRAAC))
		print("NoRAAC genes in single?")
		print(sum(NoRAAC %in% single))

		print(paste(length(intergenes_for_all_mirnas), "intergenes were found!", sep=" "))
		print(paste(length(single), "single genes were found!", sep=" "))
		print(paste(length(NoRAAC), "homoeologos genes were found!", sep=" "))

		pwf_RAAC <- c()
		pwf_NoRAAC <- c()
		pwf_single <- c()
		genes_not_found <- c()
		
		print("Selecting and filtering genes")
		i <- 0
		pb <- txtProgressBar(min=i, max=length(unique(gene2cat[,1])), style=3)
		unique_gene2cat <- unique(gene2cat[,1])
		for (item in unique_gene2cat) {
			pwf_RAAC <- c(pwf_RAAC, ifelse(item %in% intergenes_for_all_mirnas, 1, 0))
			pwf_NoRAAC <- c(pwf_NoRAAC, ifelse(item %in% NoRAAC, 1, 0))
			pwf_single <- c(pwf_single, ifelse(item %in% single, 1, 0))
			genes_not_found <- c(genes_not_found, ifelse(!(item %in% c(intergenes_for_all_mirnas, NoRAAC, single)), TRUE, FALSE))
			i <- i + 1
			setTxtProgressBar(pb, i)
		}
		print(unique_gene2cat[genes_not_found])
		
		cat("\n")
		print("Boolean vectors are unique?")
		print(sum(apply(cbind(pwf_RAAC, pwf_NoRAAC, pwf_single), 1, sum)))
		print(length(unique_gene2cat))
		print(length(unique_gene2cat) == sum(apply(cbind(pwf_RAAC, pwf_NoRAAC, pwf_single), 1, sum)))

		pwf_RAAC <- cbind(pwf_RAAC, rep(0, length(pwf_RAAC)), rep(0.5, length(pwf_RAAC)))
		rownames(pwf_RAAC) <- unique_gene2cat
		colnames(pwf_RAAC) <- c("DEgenes", "bias.data", "pwf")
		pwf_RAAC <- as.data.frame(pwf_RAAC)
		pwf_NoRAAC <- cbind(pwf_NoRAAC, rep(0, length(pwf_NoRAAC)), rep(0.5, length(pwf_NoRAAC)))
		rownames(pwf_NoRAAC) <- unique_gene2cat
		colnames(pwf_NoRAAC) <- c("DEgenes", "bias.data", "pwf")
		pwf_NoRAAC <- as.data.frame(pwf_NoRAAC)
		pwf_single <- cbind(pwf_single, rep(0, length(pwf_single)), rep(0.5, length(pwf_single)))
		rownames(pwf_single) <- unique_gene2cat
		colnames(pwf_single) <- c("DEgenes", "bias.data", "pwf")
		pwf_single <- as.data.frame(pwf_single)
		
		print("GO functional term analysis...")
		GO_RAAC <- goseq(gene2cat=gene2cat, genome=genome_file, id=unique_gene2cat, method="Hypergeometric", pwf=pwf_RAAC, use_genes_without_cat=T)
		GO_NoRAAC <- goseq(gene2cat=gene2cat, genome=genome_file, id=unique_gene2cat, method="Hypergeometric", pwf=pwf_NoRAAC, use_genes_without_cat=T)
		GO_single <- goseq(gene2cat=gene2cat, genome=genome_file, id=unique_gene2cat, method="Hypergeometric", pwf=pwf_single, use_genes_without_cat=T)
		GO_result <- SimpleList(RAAC=GO_RAAC, NoRAAC=GO_NoRAAC, single=GO_single)

		if (cmd == "plot"){
			sig_terms <- c()
			pvalues_colnames <- c(names(GO_result@listData))
			for (item in GO_result@listData) {
				sig_terms <- c(sig_terms, item[item$ontology == "BP" & item$over_represented_pvalue <= 0.01,"term"])
			}
			sig_terms <- sig_terms[!is.na(sig_terms)]
			sig_terms <- unique(sig_terms)
			pvalues <- c()
			for(term in sig_terms) {
				row <- c()
				for(item in GO_result@listData) {
					row <- c(row, item[item$ontology %in% "BP" & item$term %in% term,"over_represented_pvalue"])
				}
				pvalues <- rbind(pvalues, row)
			}
			pvalues <- -log10(pvalues)
			colnames(pvalues) <- pvalues_colnames
			rownames(pvalues) <- sig_terms
			par(mar=c(25,5,2,2))
			barplot(t(pvalues), col=c("blue", "red", "green"), beside=T, las=2, cex.names=0.8, ylab="-log10(pvalue; X >= x)")
		#	text(cex=0.7, x=bp-.25, y=-1.25, rownames(pvalues), xpb=T, srt=45)
			legend("topright", legend=c("RAAC Homeologs", "No-RAAC Homeologs", "Singletons" ), fill=c("blue", "red", "green"))
		}
		else {
			GO_result
		}
})

# Functions used to match time scaled mirnas families to its intersection quotient (QI)

setGeneric("select.species_mirna", function(mir2seeds, time_scaled_mirnas) standardGeneric("select.species_mirna"))
setMethod("select.species_mirna", signature("matrix", "data.frame"),
	function(mir2seeds, time_scaled_mirnas) {
		mirnas <- as.vector(mir2seeds[,1])
		selected_time_scaled_mirnas <- time_scaled_mirnas[time_scaled_mirnas[,1] %in% mirnas,]
		max_time_scaled_position <- max(selected_time_scaled_mirnas[,3])
		species_specific_mirnas <- mirnas[!(mirnas %in% time_scaled_mirnas[,1])]
		rows <- c()
		for(item in species_specific_mirnas){
			rows <- rbind(rows, c(item, "specific sp", max_time_scaled_position+1))
		}
		selected_time_scaled_mirnas <- rbind(selected_time_scaled_mirnas, rows)
		selected_time_scaled_mirnas
})

setGeneric("match.mirnas2QI", function(intersections, mir2seeds, colname, time_scaled_mirnas) standardGeneric("match.mirnas2QI"))
setMethod("match.mirnas2QI", signature("matrix", "matrix", "character", "data.frame"),
	function(intersections, mir2seeds, colname, time_scaled_mirnas) {
		rownames(intersections) <- gsub("\\.data\\..*", "", rownames(intersections))
		qis <- c()
		for (i in 1:length(mir2seeds[,1])) {
			qis <- c(qis, intersections[rownames(intersections) == mir2seeds[i,2], colname])
		}
		print(qis)
		if (length(mir2seeds[,1]) != length(qis)) {
			stop("mir2seeds[,1] vector length does not match Intersection quotient qis vector length")
		}
		else {
			mir2seeds <- cbind(mir2seeds, qis)
		}
		time_scaled_position <- c()
		max_time_scaled_position <- max(time_scaled_mirnas[,3]) + 1
		for (i in 1:length(mir2seeds[,1])) {
			if (mir2seeds[i,1] %in% time_scaled_mirnas[,1]) {
				time_scaled_position <- c(time_scaled_position, time_scaled_mirnas[time_scaled_mirnas[,1] == mir2seeds[i,1],3])
			}
			else {
				time_scaled_position <- c(time_scaled_position, max_time_scaled_position)
			}
		}
		print(time_scaled_position)
		if (length(mir2seeds[,1]) != length(time_scaled_position)) {
			stop("mir2seeds[,1] vector length does not match time_scaled_position vector length")
		}
		else {
			mir2seeds <- cbind(mir2seeds, time_scaled_position)
		}
		as.data.frame(mir2seeds)

})

.snippet <- function(){
genes <- genes[grepl("\\.[LSX]$", genes)]
genes
sub(".+\\.","",genes)
split(genes, sub(".+\\.","",genes))
l = split(genes, sub(".+\\.","",genes))
names(l)
class(l)
l = split(sub("\\.[LSX]$","",genes), sub(".+\\.","",genes))
lapply(l, head)
sapply(l, length)
length(intersect(l[[1]], l[[2]]))
intersect(l[[1]], l[[2]])
Reduce(l, intersect)
Reduce(intersect, l)
commonGenes <- Reduce(intersect, l)
length(commonGenes)
history()
results <- list()
history()
fileName <- "UUUCAGU.data.fu"
results[[fileName]] <- commonGenes
results
names(results)
split(names(results), sub(".+\\.","",names(results)))
}
