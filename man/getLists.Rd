% File man/getLists.Rd
\name{getLists}
\alias{getLists}
\title{Returns a list of lists of homoeologs, singletons, CRH, URH & CURH genes}
\description{
	Takes three character vectors: total genes, path where *.data.f* files are, and the stage ["fu"|"fd"]. It returns the homoeologue, singleton, CRH, URH, NRH, SIN & the genes and #miRNA families table for every one vector.

	homoeologs: Homologue genes of an hybrid specie, which come from different parental line.
	singleton: Homologue genes of an hybrid specie, which come from different parental line but one gene has lost its partner.
	CRH: Homoeologue genes which are regulated as a couple by at least one microRNA family.
	URH: Homoeologue genes which are not regulated as a couple by at least one microRNA family.
	NRH: Homoeologue genes which are not regulated by any microRNA family.
	SIN: Singleton genes which are regulated by at least one microRNA family.
	Also it returns the SIN, CRH and URH tables of gene name and # regulatory microRNAs.
	Total genes, homoeologues and singletons gene names also are returned.
}
\usage{
	getLists(total_genes = tgenes, genes = miRNA_targets)
}
\arguments{
	\item{total_genes}{A vector character of the total genes}
	\item{genes}{A character vector of the predicted miRNA taget genes}
}
\seealso{
	\code{\link{}}
}
\examples{
	lists <- getLists(total_genes = tgenes, genes = genes)
}
\keyword{file}
