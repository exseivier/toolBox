% File man/match.genes2UTRlen.Rd
\name{match.genes2UTRlen}
\alias{match.genes2UTRlen}
\title{Matches gene names of intergenes table to UTR length of genes2UTRlen table}
\description{
	Takes two matrix of data: (1) genes2UTRlen which contains the gene name and the length of its 3' UTR in bp.
	(2) intergenes is another matrix which contains the frequency of observing this gene regulated by one or
	mor miRNAs families. This dataset must contain the gene identifier as rownames and should match with the
	names of the genes in genes2UTRlen except for the .L or .S tag
}
\usage{
	match.genes2UTRlen(intergenes=intergenes, genes2UTRlen=genes2UTRlen)
}
\arguments{
	\item{intergenes}{This dataset contains the frequency of observing the homoeolog genes be regulated by at least by one miRNA as a couple}
	\item{genes2UTRlen}{This dataset contains the name of the genes and the length of its 3' UTR.}
}
\seealso{
	\code{\link{Pending}}
}
\examples{
	intergenes <- read.delim("file", header=T, row.names=1, sep="\t")
	genes2UTRlen <- read.delim("file2.txt", header=T, row.names=1, sep="\t")
	intergenes <- as.matrix(intergenes)
	genes2UTRlen <- as.matrix(genes2UTRlen)
	final <- match.genes2UTRlen(intergenes, genes2UTRlen)
}
\keyword{file}
