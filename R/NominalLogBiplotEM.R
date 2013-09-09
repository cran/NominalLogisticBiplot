# file NominalLogisticBiplot/R/NominalLogBiplotEM.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#This function is an alternated algorith for calculating en two steps the row coordinates,
# or also known as ability (E-step), and parameters of each model for the variables (M-step).
# The algorith uses multidimensional Gauss-Hermite quadrature and Ridge Regression to solve the
# separation problem in logistic regression.
#----------------------Parameters--------------
  #x: matrix with the nominal data 
  #dim: dimension of our solution or reduced space
  #penalization: value to correct the separation problem through the ridge regression
  #cte : it will be true if the model has an independent term.
  #tol : value to decide if the algorith should continue
  #maxiter : value to decide if the algorith should continue
  #Plot: Boolean parameter to decide if we want to plot the ability(rows)
  #showResults: boolean parameter if we want to see values in each iteration of the process.
NominalLogBiplotEM <- function(x, dim = 2, nnodos = 10, tol = 1e-04, maxiter = 100, penalization = 0.2,Plot=FALSE,showResults=FALSE) {
	################## Nodes and weights from the Gauss-Hermite quadrature
	Q = multiquad(nnodos, dim)
	X = Q$X
	A = Q$A
	q = dim(X)[1]
  p = dim(x)[2]
  Ncats=ColMax(x)
  Maxcat=max(Ncats)
	G = NominalMatrix2Binary(x)
  s = dim(G)[1]
	n = dim(G)[2]
	#Array to keep the models for each variable using Ridge Regression
	VariableModels = 0
	# Inicial parameters
	par = array(0, c(Maxcat-1, dim + 1, p))
	corr=afc(G,neje=dim)
	ability=corr$RowCoordinates[,1:dim]
	logLikold=0
	for (j in 1:p){
  	model = RidgeMultinomialRegression(x[, j], ability, penalization = penalization, tol = tol, maxiter = maxiter,showIter = showResults)  	
  	VariableModels=cbind(VariableModels,model)
  	par[1:(Ncats[j]-1),,j] =model$beta
  	logLikold=logLikold+model$logLik
  }
  VariableModels=VariableModels[,2:(p+1)]
  
	if (Plot){
  	plot(ability[,1],ability[,2])
    text(ability[,1],ability[,2],1:q, pos=4, offset=0.2)
  }
	#################### Inicialization of the parameters used in each iteration
	error = 1
	iter = 0
	while ((error > tol) & (iter < maxiter)) {
		# E-step
		iter = iter + 1
		PT= EvalPolylogist(X, par, Ncats)
		L = matrix(1, s, q)
		for (l in 1:s)
     for (k in 1:q)
       L[l, k] = prod(PT[k,]^G[l,])
		Pl = L %*% A
		ability = matrix(0, s, dim)
		for (l in 1:s)
     for (j in 1:dim){
			for (k in 1:q)
         ability[l, j] = ability[l, j] + X[k, j] * L[l, k] * A[k]
			ability[l, j] = ability[l, j]/Pl[l]
		}
  	###################### M-step  -  Calculation of the parameters
  	logLik=0
    for (j in 1:p){
 	  	model = RidgeMultinomialRegression(x[, j], ability, penalization = penalization, tol = tol, maxiter = maxiter,showIter = showResults)  	
 	  	VariableModels[,j] = model
    	par[1:(Ncats[j]-1),,j] =model$beta
    	logLik=logLik+model$logLik
    }
		error = abs((logLik - logLikold)/logLik)
		if(showResults) print(c(iter, logLik, error))
		logLikold = logLik
		if(showResults) print(paste("iter:",iter," error:",error,sep=""))
	}
	if (Plot)
    	for (j in 1:p){
      	dev.new()
      	plot(ability[,1],ability[,2],col=x[,j])
        text(ability[,1],ability[,2],1:q,col=x[,j])
      }
	d = sqrt(rowSums(par^2) + 1)
	model = list()
	model$RowCoordinates = ability
	model$ColumnModels = VariableModels
	class(model) = "nominal.logistic.biplot.EM"
	return(model)
}
