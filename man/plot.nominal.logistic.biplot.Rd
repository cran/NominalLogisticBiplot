\name{plot.nominal.logistic.biplot}
\alias{plot.nominal.logistic.biplot}
\title{
Graphical representation of a Nominal Logistic Biplot.
}
\description{
Plotting a Nominal Logistic Biplot. There are parameters related to the way in which the biplot is plotted. All the posible parameters have default values
}
\usage{
## S3 method for class 'nominal.logistic.biplot'
\method{plot}{nominal.logistic.biplot}(x, planex = 1, planey = 2,
QuitNotPredicted = TRUE, ReestimateInFocusPlane = TRUE,
proofMode = FALSE, AtLeastR2 = 0.01, xlimi = -1.5, xlimu = 1.5,
ylimi = -1.5, ylimu = 1.5, linesVoronoi = FALSE, ShowAxis = TRUE,
PlotVars = TRUE, PlotInd = TRUE, LabelVar = TRUE, LabelInd = TRUE,
CexInd = NULL, CexVar = NULL, ColorInd = NULL, ColorVar = NULL,
SmartLabels = FALSE, PchInd = NULL, PchVar = NULL,
LabelValuesVar = NULL, ShowResults = FALSE,\dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of the class nominal.logistic.biplot.
}
  \item{planex}{
Dimension for X axis.
}
  \item{planey}{
Dimension for Y axis.
}
  \item{QuitNotPredicted}{
  Should the non-predicted categories  be represented on the graph?
}
  \item{ReestimateInFocusPlane}{
  Shuld the item parameters be reestimated using only the dimensiona of the plot.? If FALSE the values of the parameters for other dimensions are set to 0. Default is FALSE}
  \item{proofMode}{
 Should each variable be plotted on a separate plot? If FALSE, a single plot with a legend for identifying each variable is made.
}
  \item{AtLeastR2}{
It establishes the cutting value to plot a variable attending to its Nagelkerke R^2  value. A variable is plotted if its R^2 is higher than this value.
}
  \item{xlimi}{
Minimum value on the x-axis.
}
  \item{xlimu}{
Maximum value on the x-axis.
}
  \item{ylimi}{
Minimum value on the y-axis.
}
  \item{ylimu}{
Maximum value on the y-axis.
}
  \item{linesVoronoi}{
Should the tesselation be plotted.? Default is FALSE and only the category points are plotted for a better reading of the plot. 
}
  \item{ShowAxis}{
  Should the axis be shown?
}
  \item{PlotVars}{
  Should the variables (items) be pplotted?
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{LabelVar}{
Should the variable labels be shown?
}
  \item{LabelInd}{
Should the individual labels be shown?
}
  \item{CexInd}{
Size of the individual points. It can be an array with the cex information for each row.
}
  \item{CexVar}{
Size of the category points. It can be an array with the cex information for each variable.
}
  \item{ColorInd}{
Color of the individual points. It can be an array with the color information for each row.
}
  \item{ColorVar}{
Color for the variables. It can be an array with the color information for each variable.
}
  \item{SmartLabels}{
Should the text labels be printed accordind to its position on the plot?.            
}
  \item{PchInd}{
Symbol for the individuals. It can be an array with the pch information for each row.
}
  \item{PchVar}{
Symbol for the variables. It could be an array with the pch information for each variable.
}
  \item{LabelValuesVar}{
List with the text labels for all the variables. If NULL, initial labels are used.
}
  \item{ShowResults}{
Should the results of the proccess of calculating the prediction regions be shown?
}
  \item{\dots}{
Additional parameters to plot.
}
}
\details{
The function without parameters plots the nominal.logistic.biplot object with labels in the original data and default values for colors, symbols and sizes for points and lines. Other values of colors, symbols and sizes can be supplied. A single value applies to all the points but an array with different values can be used to improve the undestanding of the plot.-
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
plot(nlbo,QuitNotPredicted=TRUE,ReestimateInFocusPlane=TRUE,
planex = 1,planey = 2,proofMode=TRUE,LabelInd=TRUE,
AtLeastR2 = 0.01,xlimi=-1.5,xlimu=1.5,ylimi=-1.5,
ylimu=1.5,linesVoronoi = TRUE,SmartLabels = FALSE,
PlotInd=TRUE,CexInd = c(0.6,0.7,0.5,0.4,0.5,0.6,0.7)
,PchInd = c(1,2,3,4,5,6,7),ColorInd="black",PlotVars=TRUE,
LabelVar = TRUE,PchVar = c(1,2,3,4,5),
ColorVar = c("red","black","yellow","blue","green")
,ShowResults=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

