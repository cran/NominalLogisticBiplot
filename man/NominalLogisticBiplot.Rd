\name{NominalLogisticBiplot}
\alias{NominalLogisticBiplot}
\title{Nominal Logistic Biplot for polytomous data}
\description{
Function that calculates the parameters of the Nominal Logistic Biplot according to Hernandez-Sanchez & Vicente-Villardon (2013). 
}

\usage{
NominalLogisticBiplot(datanom, sFormula = NULL, numFactors = 2,
coordinates = "EM", rotation = "varimax", metfsco = "EAP",
nnodos = 10, tol = 1e-04, maxiter = 100, penalization = 0.1,
cte = TRUE, show = FALSE)
}
\arguments{
  \item{datanom}{
The data set, it can be a \emph{matrix} with integers or a \emph{data frame} with factors. All variables have to be nominal. 
}
  \item{sFormula}{
This parameter follows the unifying interface for selecting variables from a data frame for a plot, test or model. The most common formula es of type y ~ x1+x2+x3. It has a default   value of NULL if not specified.
}
  \item{numFactors}{
Number of dimensions of the solution. It should be lower than the number of variables. It has a default value of 2.
}
  \item{coordinates}{
This parameter can be: "EM", "ACM", "MIRT" or "PCOA". Method to compute the row coordinates.
}
  \item{rotation}{
Rotation method to used with "MIRT" option in  "coordinates". No effect fror other options.
}
  \item{metfsco}{
Calculation method for the fscores with  "MIRT" option in  "coordinates". No effect fror other options.
}
  \item{nnodos}{
Number of nodes for gauss quadrature in the EM algorithm.
}
  \item{tol}{
Tolerance for the EM algorithm.
}
  \item{maxiter}{
Maximum number of iterations in the EM algorithm.
}
  \item{penalization}{
Penalization for the ridge regression for each variable.
}
  \item{cte}{
Include constant in the logistic regression model. Default is TRUE.
}
  \item{show}{
Show intermediate copmputations. Default is TRUE.
}
}
\details{
The general algorithm used is essentially an alternated procedure in which parameters for rows and columns are computed in alternated steps repeated until convergence. Parameters for the rows are calculated by expectation (E-step) or by a external procedure (Multiple Correspondence Analysis or Principal Coordinates Analysis) and parameters for the columns are computed by maximization (M-step), i. e., by Nominal Logistic Regression. When the procedure for Row scores is  external, only one iteration is performed and the procedure is called "External Nominal Logistic Biplot". Because the aim of the biplot is the representation 


There are several options for the computation:

1.- Using the package \pkg{mirt} to obtain the row scores, i. e. using a solution obtained from a latent trait model. The column (item) parameters should be directly used by our biplot procedure but, because of the characteristics of the package that performs a default rotation after parameter estimation, we have to reestimate the item parametes to be coherent to the scores.

2.- Using our implementation of the EM algorithm alternating expected a porteriori scores and Ridge Nominal Logistic Regression for each variable.

3.- Using external coordinates for the rows taken from Multiple Correspondence Analysis or Principal Coordinates Analysis and fitting the response surfaces in just one step.


Equations that define a set of probability response surfaces (one for each category and each variable) are no longer sigmoid as in the binary case (Vicente-Villardon et al. (2006)). This means that the level curves are no longer straight lines and then, prediction of probabilities is not made by projection as in the usual linear biplots. For each variable, define a set of convex polygons that can be interpreted as "prediction" regions in the same way as in Gower & Hand (1996).
Each pair of response surfaces defined by intersect in a straight line that, projected onto the space of predictors, is the set of points in which the probability of both categories is the same. Those lines are the candidates to be the edges of the convex polygons defining the prediction regions.
}
\value{
An object of class \code{"nominal.logistic.biplot"}.  This has some components:

  \item{dataSet }{Data set of study with all the information about the name of the levels and names of the variables and individuals}
  \item{ItemsCoords }{Coordinates for the rows in the reduced space}
  \item{VariableModels }{Information of the regression resuls for each variable.}
  \item{NumFactors  }{Number of dimensions selected for the study}
  \item{Coordinates }{Type of coordinates to calculate the row positions}
  \item{Rotation }{Type of rotation if we have chosen mirt coordinates}
  \item{Methodfscores }{Method of calculation of the fscores in mirt process}
  \item{NumNodos }{Number of nodes for the gauss quadrature in EM algorith}
  \item{tol }{Cut point to stop the EM-algorith}
  \item{maxiter }{Maximum number of iterations in the EM-algorith}
  \item{penalization }{Value for the correction of the ridge regression}
  \item{cte }{Boolean value to choose if the model for each variable will have
              independent term}
  \item{show }{Boolean value to indicate if we want to see the results of our analysis}
}
\references{

Hernandez, J. C. & Vicente-Villardon, J. L. (2013) Logistic Biplots for Nominal Data. Submitted. Preprint available at : 
https://www.researchgate.net/publication/256428288_Logistic_Biplot_for_Nominal_Data?ev=prf_pub

Vicente-Villardon, J., Galindo, M.P & Blazquez-Zaballos, A. (2006), \emph{Logistic biplots},Multiple Correspondence Analysis and related methods pp. 491--509.

Demey, J., Vicente-Villardon, J. L., Galindo, M.P.  & Zambrano, A. (2008) \emph{Identifying Molecular Markers Associated With Classification Of Genotypes Using External Logistic Biplots}. Bioinformatics, 24(24), 2832-2838.

Baker, F.B. (1992): Item Response Theory. Parameter Estimation Techniques. Marcel Dekker. New York.

Gabriel, K. (1971), \emph{The biplot graphic display of matrices with application to
principal component analysis.}, Biometrika 58(3), 453--467.

Gabriel, K. R. (1998), \emph{Generalised bilinear regression}, Biometrika 85(3),
689--700.

Gabriel, K. R. & Zamir, S. (1979), \emph{Lower rank approximation of matrices by least
squares with any choice of weights}, Technometrics 21(4), 489--498.

Gower, J. & Hand, D. (1996), \emph{Biplots, Monographs on statistics and applied probability}. 54. London: Chapman and Hall., 277 pp.

Chalmers,R,P (2012). \emph{mirt: A Multidimensional Item Response Theory Package for the R Environment}. Journal of Statistical Software, 48(6), 1-29. URL http://www.jstatsoft.org/v48/i06/.
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}

\seealso{
\code{\link{NominalLogBiplotEM}}, \code{\link{afc}}, \code{\link{PCoA}}
}
\examples{

data(HairColor)
nlbo = NominalLogisticBiplot(HairColor,sFormula=NULL,
numFactors=2,coordinates="EM",penalization=0.2,show=FALSE)
nlbo

}
\keyword{biplot}
\keyword{logistic}
\keyword{models}
