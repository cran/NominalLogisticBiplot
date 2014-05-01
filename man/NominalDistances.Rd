\name{NominalDistances}
\alias{NominalDistances}
\title{
Distances between indidividuals calculated from nominal variables.
}
\description{
This function calculates the hamming distances (or similarities) among individuals from a nominal data matrix. 
}
\usage{
NominalDistances(x, similarities = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
This parameter is a matrix with the nominal variables
}
  \item{similarities}{
Boolean parameter to specify if the user wants a distances matrix or a similarities matrix. By default this parameter is FALSE, so the function calculates the distances.
}
}
\value{
The function returns a matrix with the distances or similarities
}
\references{
Boriah, S., Chandola, V. & Kumar,V.(2008) \emph{Similarity measures for categorical data: A comparative evaluation}. In proceedings of the eight SIAM International Conference on Data Mining, pp 243--254
}
\author{
Jose Luis Vicente-Villardon,Julio Cesar Hernandez Sanchez

Maintainer: Jose Luis Vicente-Villardon <villardon@usal.es>
}
\examples{

  data(HairColor)
  NominalDistances(data.matrix(HairColor))
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nominal}
\keyword{distances}% __ONLY ONE__ keyword per line
