\name{afc}
\alias{afc}
\title{
Simple Correspondence Analysis
}
\description{
This function calculates simple correspondence analysis for a data matrix. The function is used to calculate multiple correspondence analysis as a starting solution of the EM procedure.
}
\usage{
afc(x, neje = 2, alfa = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A frequency matrix or a binary matrix obtained from the original data set of nominal variables.
}
  \item{neje}{
Number of dimensions for the solution
}
  \item{alfa}{
Biplot weight for rows and columns. 1 means rows in principal coordinates and columns in standard coordinates, 0 means rows in standard coordinates and columns in principal coordinates.
}
}

\value{
An object of class \code{"afc.sol"}.This has some components:

  \item{RowCoordinates }{Coordinates for the rows in the reduced dimension space}
  \item{ColCoordinates }{Coordinates for the variables in the reduced dimension space}
  \item{RowContributions }{Contributions of the dimensions to explain the inertia of each row}
  \item{ColContributions }{Contributions of the dimensions to explain the inertia of each row}
  \item{Inertia }{Inertia for each dimension}
  \item{Eigenvalues }{Eigenvalues}
}
\references{
BENZECRI, J.P. (1973) \emph{L'analyse des Donnees}. Vol. 2. L'analyse des correspondences. Dunod. Paris.

}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{NominalMatrix2Binary}}
}
\examples{
  
  data(HairColor)
  G = NominalMatrix2Binary(data.matrix(HairColor))
	mca=afc(G,neje=2)
  mca
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{correspondence}
\keyword{analysis}% __ONLY ONE__ keyword per line
