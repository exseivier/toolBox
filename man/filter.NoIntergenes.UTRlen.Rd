% File man/filter.NoIntergenes.UTRlen.Rd
\name{filter.NoIntergenes.UTRlen}
\alias{filter.NoIntergenes.UTRlen}
\title{Picks up the No Intergenes (no-RAAC genes) from the whole set of genes}
\description{
	Takes a matrix vector which contains a numeric column (the frequency of observing this pair of homoeologs being regulated by one or more miRNA families), and a data frame which contain the names of genes matching its 3'UTR length and returns a matrix containing the names of the genes, a column with the 3' UTR length, and the corresponding number of the first column from the input matrix.
}
\usage{
	filter.NoIntergenes.UTRlen(intergenes=intergenes, genes2UTRlen=genes2UTRlen)
}
\arguments{
	\item{intergenes}{The matrix with the gene names and the frequency observed of being regulated by one or more miRNA families as a couple}
	\item{genes2UTRlen}{The data frame with the names of genes matching its 3' UTR length}
}
\seealso{
	\code{\link{match.genes2UTRlen}}
}
\examples{
	NoRAAC <- filter.NoIntergenes.UTRlen(intergenes=matrix, genes2UTRlen=data.frame)
}
\keyword{file}
