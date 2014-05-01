\name{Nominal2Binary}
\alias{Nominal2Binary}
\title{
Transformation of a nominal variable into a binary indicator matrix
}
\description{
This function transforms a nominal variable into a binary matrix with as many
  colums as categories. Each row of the matrix has
   a value of 1 for the corresponding level of the category and 0 elsewhere.
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

  \item{Z }{The binary indicator matrix asociated to the nominal variable}

}

\author{
Jose Luis Vicente-Villardon,Julio Cesar Hernandez Sanchez

Maintainer: Jose Luis Vicente-Villardon <villardon@usal.es>
}
\examples{
  
  data(HairColor)
  Nominal2Binary(as.numeric(HairColor[,1]))
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nominal}
\keyword{binary}% __ONLY ONE__ keyword per line
