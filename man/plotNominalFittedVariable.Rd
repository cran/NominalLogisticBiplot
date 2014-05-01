\name{plotNominalFittedVariable}
\alias{plotNominalFittedVariable}
\title{
Function for plotting in the reduced space an unordered and fitted categorical variable.
}
\description{
Graphical representation of a polytomous unordered variable previously fitted in the reduced space,
according to the Nominal Logistic Biplot theory. It can be choosen some parameters related to the
way in which the variable is plotted.
}
\usage{
plotNominalFittedVariable(nameVar, numcateg, beta, varstudyC, rowCoords,
levelsVar = NULL, numFactors = 2, planex = 1, planey = 2, xi = -3.5, xu = 3.5,
yi = -3.5, yu = 3.5, CexVar = 0.7,ColorVar = "blue", PchVar = 0.7,
addToPlot = FALSE, QuitNotPredicted = TRUE, ShowResults = TRUE,
linesVoronoi = TRUE, LabelVar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nameVar}{
Name of the variable to be plotted.
}
  \item{numcateg}{
Number of categories of the variable.
}
  \item{beta}{
Estimated coefficients matrix.
}
  \item{varstudyC}{
Values of the categorical variable to be plotted. It should be a factor with information
about a nominal variable, i.e., an unordered variable.}
  \item{rowCoords}{
Estimation coordinates for the individuals in the spanned space.
}
  \item{levelsVar}{
Vector with the labels for each level of the variable.
}
  \item{numFactors}{
Dimension of the reduced space.
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
Should the tesselation be plotted.? Default is FALSE and only the category
points are plotted for a better reading of the plot
}
  \item{LabelVar}{
Should the variable labels be shown?
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
data(Env)
nlbo = NominalLogisticBiplot(Env,sFormula=NULL,
numFactors=2,method="EM",penalization=0.2,show=FALSE)
nameVar = nlbo$dataSet$ColumNames[1]
numcateg = 4
beta = nlbo$VariableModels[,1]$beta
Nagelkerke = nlbo$VariableModels[,1]$Nagelkerke 
varstudyC = as.matrix(as.numeric(Env[,1]))
rowCoords = nlbo$RowsCoords
levelsVar = c("M1","M2","M4","M5")
plotNominalFittedVariable(nameVar,numcateg,beta,varstudyC,rowCoords,levelsVar=NULL,
            numFactors=2,planex = 1,planey = 2,xi=-3.5,xu=3.5,yi=-3.5,yu=3.5,
            CexVar=0.7,ColorVar="blue",PchVar=0.7,addToPlot=FALSE,
            QuitNotPredicted=TRUE,ShowResults=TRUE,linesVoronoi=TRUE,LabelVar=TRUE)

}
\keyword{nominal }
\keyword{fitting }
\keyword{plot }
