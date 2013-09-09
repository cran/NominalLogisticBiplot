\name{hermquad}
\alias{hermquad}
\title{
Gauss-Hermite quadrature
}
\description{
Computes the Hermite Quadrature weights for a set of grid points
}
\usage{
hermquad(N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Number of nodes
}
}

\value{
An object of class \code{"GaussQuadrature"}.  This has the components:

  \item{X }{Coordinates of the nodes}
  \item{W }{Weights asociated to each node}

}
\references{
Stroud, A.H. and Secrest, D. (1966) \emph{Gaussian Quadrature Formulas}, Englewood Cliffs, NJ: Prentice-Hall.

 Hildebrand,F. B. (1987) \emph{Intoduction to Numerical Analysis 2nd Ed}, Dover Publications, New York,  page 385
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}

\examples{
  
  hermquad(10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{gauss}
\keyword{hermite}
\keyword{quadrature}
