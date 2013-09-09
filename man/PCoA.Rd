\name{PCoA}
\alias{PCoA}
\title{
Principal Coordinates Analysis
}
\description{
This function calculates principal coordinates analysis using a distante
  matrix among a set of objets.
}
\usage{
PCoA(dis, r = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dis}{
Distance matrix between a set ob objects.
}
  \item{r}{
Number of dimensions for the solution.
}
}

\value{
An object with has some components:

  \item{EigenValues }{Eigenvalues of the inner products matrix}
  \item{Inertia }{Variance (Inertia) accounted for each dimension}
  \item{RowCoordinates }{Coordinates for the rows in the reduced space}
  \item{RowQualities }{Qualities of representation of the objects. Squared cosines between the points (vectors) in the full space and the points in the reduced space. Values near 1 indicate good quality}
}
\references{
Gower,J.C. (1966) \emph{Some distance properties of latent root and vector methods used in multivariate analysis}. Biometrika, 53, 325--338.
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}

\seealso{
\code{\link{NominalDistances}}
}
\examples{
  
  data(HairColor)
  dis = NominalDistances(data.matrix(HairColor))
  PCoA(dis,2)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PCoA}

