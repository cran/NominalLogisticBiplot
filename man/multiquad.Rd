\name{multiquad}
\alias{multiquad}
\title{
Multidimensional Gauss-Hermite quadrature
}
\description{
This function computes the gauss-hermite quadrature in more than
  one dimension. 
}
\usage{
multiquad(nnodos, dims)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nnodos}{
Number of nodes.
}
  \item{dims}{
Number of dimensions of the quadrature
}
}

\value{
An object of class \code{"MultiGaussQuadrature"}.  This has the components:

  \item{X }{Coordinates of the nodes}
  \item{A }{Weights asociated to each node}
}
\references{
Jackel, P. (2005) \emph{A note on multivariate Gauss-Hermite quadrature}

http://www.pjaeckel.webspace.virginmedia.com/ANoteOnMultivariateGaussHermiteQuadrature.pdf
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}

\seealso{
\code{\link{hermquad}}
}
\examples{
  
  multiquad(10,2)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{gauss}
\keyword{quadrature}% __ONLY ONE__ keyword per line
