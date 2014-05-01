\name{NominalLogisticBiplot-package}
\alias{NominalLogisticBiplot-package}
\docType{package}
\title{
Nominal Logistic Biplot representations for polytomous data
}
\description{
Analysis of a matrix of polytomous items using Nominal Logistic Biplots (NLB) according to Hernandez-Sanchez & Vicente-Villardon (2013). The NLB procedure extends the binary logistic biplot to nominal (polytomous) data. The  individuals are represented as points on a plane and the  variables are represented as convex prediction regions rather than vectors as in a classical or binaly biplot. Using the methods from the
  Computational Geometry, the set of prediction regions is converted to a set of points
  in such a way that the prediction for each individual is established by its closest
  "category point". Then interpretation is based on distances rather than on projections.
  In this package we implement the geometry of such a representation and construct computational algorithms
  for the estimation of parameters and the calculation of prediction regions}
\details{
\tabular{ll}{
Package: \tab NominalLogisticBiplot\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-08-05\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon
Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\keyword{package}
\seealso{
\code{\link{NominalLogisticBiplot}},\code{\link{NominalLogBiplotEM}},\code{\link{multiquad}},\code{\link{summary.nominal.logistic.biplot}},\code{\link{plot.nominal.logistic.biplot}}
}
\examples{
data(HairColor)
nlbo = NominalLogisticBiplot(HairColor,sFormula=NULL,numFactors=2,
method="EM",penalization=0.2,show=FALSE)
summary(nlbo)
plot(nlbo,QuitNotPredicted=TRUE,ReestimateInFocusPlane=TRUE,
planex = 1,planey = 2,proofMode=TRUE,LabelInd=TRUE,AtLeastR2 = 0.01
,xlimi=-1.5,xlimu=1.5,ylimi=-1.5,ylimu=1.5,linesVoronoi = TRUE
,SmartLabels = FALSE, PlotInd=TRUE,CexInd = c(0.6,0.7,0.5,0.4,0.5,0.6,0.7)
,PchInd = c(1,2,3,4,5,6,7),ColorInd="black",PlotVars=TRUE,LabelVar = TRUE
,PchVar = c(1,2,3,4,5),ColorVar = c("red","black","yellow","blue","green")
,ShowResults=TRUE)

}
