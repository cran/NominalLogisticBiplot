\name{summary.nominal.logistic.biplot}
\alias{summary.nominal.logistic.biplot}
\title{
Summary Method Function for Objects of Class 'nominal.logistic.biplot'
}
\description{
This function shows a summary of the principal results for the estimation
  for individuals and variables, like some Pseudo R-squared indices, the correct
  classification percentage of each regression, the logLikelihood and
  "Estimate coefficients", "Std. Error", "z value" or "Pr(>|z|)" values.
}
\usage{
## S3 method for class 'nominal.logistic.biplot'
\method{summary}{nominal.logistic.biplot}(object,completeEstim,coorInd,\dots)
}
\arguments{
  \item{object}{
This parameter keeps the nominal logistic biplot object.
}
  \item{completeEstim}{
Boolean parameter to choose if the estimated coefficients will be printed on screen.Default value is FALSE.
}
  \item{coorInd}{
Boolean parameter to choose if the individual coordinates will be printed on screen.Default value is FALSE.
}
  \item{\dots}{
Additional parameters to summary.
}
}
\details{
This function is a method for the generic function summary() for class "nominal.logistic.biplot". It can be invoked by calling summary(x) for an object x of the appropriate class.
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{NominalLogisticBiplot}}
}
\examples{
data(HairColor)
nlbo = NominalLogisticBiplot(HairColor,sFormula=NULL,
numFactors=2,method="EM",penalization=0.2,show=FALSE)
summary(nlbo)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{summary}

