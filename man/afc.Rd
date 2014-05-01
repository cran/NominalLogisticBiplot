\name{afc}
\alias{afc}
\title{
Simple Correspondence Analysis
}
\description{
This function calculates simple correspondence analysis for a data matrix.
}
\usage{
afc(x, dim = 2, alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A frequency matrix or a binary matrix obtained from the original data set of nominal variables.
}
  \item{dim}{
Number of dimensions for the solution
}
  \item{alpha}{
Biplot weight for rows and columns. 1 means rows in principal coordinates and columns in standard coordinates, 0 means rows in standard coordinates and columns in principal coordinates.
}
}

\value{
An object of class \code{"afc.sol"}.This has some components:

  \item{Title }{Title of the statistical technique}
  \item{Non_Scaled_Data }{Original data}
  \item{Minima }{vector with the minimum values for each column of the initial data matrix}
  \item{Maxima }{vector with the maximum values for each column of the initial data matrix}
  \item{Initial_Transformation }{Name of the transformation for the data}
  \item{Scaled_Data }{Scaled data according to the transformation}
  \item{nrows }{Number of rows of the data set}
  \item{ncols }{Number of columns of the data set}
  \item{dim }{Number of dimensions for the solution}
  \item{CumInertia }{Acumulated Inertia}
  \item{Scale_Factor }{Scale factor for the transformation}
  \item{RowCoordinates }{Coordinates for the individuals in the reduced dimension space}
  \item{ColCoordinates }{Coordinates for the variables in the reduced dimension space}
  \item{RowContributions }{Contributions of the dimensions to explain the inertia of each row}
  \item{ColContributions }{Contributions of the dimensions to explain the inertia of each column}
  \item{Inertia }{Inertia for each dimension}
  \item{Eigenvalues }{Eigenvalues}
}
\references{
BENZECRI, J.P. (1973) \emph{L'analyse des Donnees}. Vol. 2. L'analyse des correspondences. Dunod. Paris.

}
\author{
 Jose Luis Vicente-Villardon,Julio Cesar Hernandez Sanchez

Maintainer: Jose Luis Vicente-Villardon <villardon@usal.es>
}
\seealso{
\code{\link{NominalMatrix2Binary}}
}
\examples{
  
  data(HairColor)
  G = NominalMatrix2Binary(data.matrix(HairColor))
	mca=afc(G,dim=2)
  mca
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{correspondence}
\keyword{analysis}% __ONLY ONE__ keyword per line
