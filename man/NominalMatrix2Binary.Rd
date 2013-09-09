\name{NominalMatrix2Binary}
\alias{NominalMatrix2Binary}
\title{
Indicator matrix of a set of nominal variables.
}
\description{
Constructs the indicator matrix for a nominal variables matrix.
}
\usage{
NominalMatrix2Binary(Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A matrix with nominal variables measured for a set of individuals. Input must be a matrix with integer values.
}
}
\value{
\item{G }{The binary indicator matrix asociated to the nominal matrix}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{Nominal2Binary}}
}
\examples{

  data(HairColor)
  NominalMatrix2Binary(data.matrix(HairColor))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nominal}
\keyword{binary}% __ONLY ONE__ keyword per line
