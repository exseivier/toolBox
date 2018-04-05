% File man/getLists.Rd
\name{getLists}
\alias{getLists}
\title{Returns a list of lists of homoeologs, singletons, CRH, URH & CURH genes}
\description{
	Takes two character vector total genes and predicted target genes and returns the homoeolog, singleton, CRH, URH, CURH & the genes and #miRNA families table.

	homoeologs: Homologue genes of an hybrid specie, which come from different parental line.
	singleton: Homologue genes of an hybrid specie, which come from different parental line but one gene has lost its partner.
	CRH: Homoeologue genes which are regulated as a couple by at least one microRNA family.
	URH: Homoeologue genes which are not regulated as a couple by at least one microRNA family.
	CURH: A mix of CRH & URH
	genes and #miRNA families table: It is a table composed by 1 column which contains the number of microRNAs that regulates a determinated gene.
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
