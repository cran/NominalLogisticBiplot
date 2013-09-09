\name{Nominal2Binary}
\alias{Nominal2Binary}
\title{
Transformation from a nominal variable into a binary matrix
}
\description{
This function transfomates a nominal variable into a binary matrix with as many
  colums as categories. Each row of the matrix has
  only a single one for the corresponding level of the category.
}
\usage{
Nominal2Binary(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A vector containing the values of nominal variable measured on a set of individuals- The values must be integers starting at 1.
}
}
\value{
An object of type matrix:

  \item{Z }{The binary matrix asociated to the nominal variable}

}

\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\examples{
  
  data(HairColor)
  Nominal2Binary(as.numeric(HairColor[,1]))
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nominal}
\keyword{binary}% __ONLY ONE__ keyword per line
