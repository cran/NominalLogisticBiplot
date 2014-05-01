\name{NominalLogBiplotEM}
\alias{NominalLogBiplotEM}
\title{
Alternated EM algorithm for Nominal Logistic Biplots
}
\description{
This function computes, with an alternated algorithm, the row and column parameters of a Nominal Logistic Biplot for polytomous data. The row coordinates (E-step) are computed using multidimensional Gauss-Hermite quadratures and Expected \emph{a posteriori} (EAP) scores and parameters for each variable or items (M-step)using Ridge Nominal Logistic Regression to solve the separation problem present when the points for different categories of a variable are completely separataed on the representation plane and the usual fitting methods do not converge. The separation problem is present in almost avery data set for which the goodness of fit is high. 
} 
\usage{
NominalLogBiplotEM(x, dim = 2, nnodos = 10, tol = 1e-04,
maxiter = 100, penalization = 0.2,initial=1,alfa=1, Plot = FALSE,
showResults = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix with the nominal data. The matrix must be in numerical form. 
}
  \item{dim}{
Dimension of the solution
}
  \item{nnodos}{
Number of nodes for the multidimensional Gauss-Hermite quadrature
}
 \item{tol}{
Value to stop the process of iterations.
}
  \item{maxiter}{
Maximum number of iterations in the process of solving the regression coefficients.
}
  \item{penalization}{
Penalization used in the diagonal matrix to avoid singularities.
}
  \item{initial}{
Value to decide the method(1-Correspondence analysis, 2-Mirt) that calculates the initial abilities values for the individuals. 
}
  \item{alfa}{  
If initial parameter method is correspondence analysis, this parameter determines the weight for rows and columns.
}
  \item{Plot}{
Boolean parameter to plot the row coordinates
}
  \item{showResults}{
Boolean parameter to show all the information about the iterations.
}
}

\value{
An object of class \code{"nominal.logistic.biplot.EM"}.This has components:

  \item{RowCoordinates }{Coordinates for the individuals in the reduced space}
  \item{ColumnModels }{List with information about the Nominal Logistic Models calculated for each variable including: estimated parameters with covariances and standard errors, log-likelihood, deviances, percents of correct classifications, pvalues and pseudo-Rsquared measures}

}
\references{
  Bock,R. & Aitkin,M. (1981),\emph{Marginal maximum likelihood estimation of item parameters: Aplication of an EM algorithm}, Phychometrika 46(4), 443-459.

Gabriel, K. R. (1998). Generalised bilinear regression. Biometrika, 85(3), 689-700.
  
  Vicente-Villardón, J. L., Galindo Villardón, M. P., & Blázquez Zaballos, A. (2006). Logistic biplots. Multiple correspondence analysis and related methods. London: Chapman & Hall, 503-521.
  
  Gabriel, K. R., & Zamir, S. (1979). Lower rank approximation of matrices by least squares with any choice of weights. Technometrics, 21(4), 489-498.
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{polylogist}}, \code{\link{multiquad}}
}
\examples{
    data(HairColor)
    data = data.matrix(HairColor)
    xEM = NominalLogBiplotEM(data, dim = 2,showResults = FALSE)
    xEM
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{EM}
\keyword{algorithm}
