\name{plotNominalVariable}
\alias{plotNominalVariable}
\title{
Function for plotting in the reduced space an unordered categorical variable.
}
\description{
Graphical representation of a polytomous unordered variable in the reduced space,
according to the Nominal Logistic Biplot theory. Inside the function, the estimations needed for the variable
will be done. It can be choosen some parameters related to the way in which the variable is plotted.
}
\usage{
plotNominalVariable(nameVar, nominalVar, estimRows, planex = 1, planey = 2,
xi = -3.5, xu = 3.5, yi = -3.5, yu = 3.5, CexVar = 0.7, ColorVar = "blue",
PchVar = 0.7, addToPlot = FALSE, QuitNotPredicted = TRUE, ShowResults = FALSE,
linesVoronoi = TRUE, LabelVar = TRUE, tol = 1e-04, maxiter = 100,
penalization = 0.1, showIter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nameVar}{
Name of the variable to be plotted.
}
  \item{nominalVar}{
Values of the categorical variable to be plotted. It should be a factor with information about a nominal variable, i.e., a variable without ordered values.
}
  \item{estimRows}{
Estimation coordinates for the individuals in the spanned space.
}
  \item{planex}{
Dimension for X axis.
}
  \item{planey}{
Dimension for Y axis.
}
  \item{xi}{
Minimum value on the x-axis.
}
  \item{xu}{
Maximum value on the x-axis.
}
  \item{yi}{
Minimum value on the y-axis.
}
  \item{yu}{
Maximum value on the y-axis.
}
  \item{CexVar}{
Size of the category points.
}
  \item{ColorVar}{
Color for the variable.
}
  \item{PchVar}{
Symbol for the variable.
}
  \item{addToPlot}{
  Should the graph be added to an existing representation?
}
  \item{QuitNotPredicted}{
  Should the non-predicted categories  be represented on the graph?
}
  \item{ShowResults}{
Should the results of the proccess of calculating the prediction regions be shown?
}
  \item{linesVoronoi}{
Should the tesselation be plotted.? Default is FALSE and only the category points are plotted for a better reading of the plot
}
  \item{LabelVar}{
Should the variable labels be shown?
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
  \item{showIter}{
Boolean parameter to show the information about the iterations.
}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{polylogist}}
}
\examples{
data(HairColor)
nlbo = NominalLogisticBiplot(HairColor,sFormula=NULL,
numFactors=2,method="EM",penalization=0.2,show=FALSE)
nameVar = nlbo$dataSet$ColumNames[2]
nominalVar = HairColor[,2]
estimRows = nlbo$RowsCoords 
plotNominalVariable(nameVar,nominalVar,estimRows,planex = 1,planey = 2,
    xi=-1.5,xu=1.5,yi=-1.5,yu=1.5,CexVar=0.7,ColorVar="blue",PchVar=0.7,
    addToPlot=FALSE,QuitNotPredicted=TRUE,ShowResults=TRUE,
    linesVoronoi=TRUE,LabelVar=TRUE,tol = 1e-04, maxiter = 100,
    penalization = 0.3,showIter = FALSE)

}
\keyword{nominal }
\keyword{plot }
