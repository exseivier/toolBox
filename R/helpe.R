
eliminateHeavyDuplicates <- function(names, times) {
	tmp_names <- sort(names)
	substracting_names <- sort(names)
	for(i in 1:times) {
		substracting_names <- substracting_names[duplicated(substracting_names)]
	}
	tmp_names <- tmp_names[!tmp_names %in% substracting_names]
	tmp_names
}
