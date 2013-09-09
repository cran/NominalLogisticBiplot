# file NominalLogisticBiplot/R/afc.R
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

#Function that calculates simple correspondence analysis for a data set.
#----------------------Parameters--------------
  #x: data set with nominal variables
  #neje: number of dimensions for the solution
  #alfa: this parameter is used to calculate row and column markers
afc <- function(x,neje=2, alfa=1){
    n = dim(x)[1]
    p = dim(x)[2]
    nt = sum(x)
    x = x/nt
    dr = matrix(rowSums(x), n, 1)
    dc = matrix(colSums(x), p, 1)
    x = x - dr %*% t(dc)
    Dr = diag(1, n, n)
    diag(Dr) <- 1/sqrt(dr)
    Dc = diag(1, p, p)
    diag(Dc) <- 1/sqrt(dc)
    x = Dr %*% x %*% Dc
    UDV = svd(x)
    r = min(c(n, p))
    d = UDV$d[1:r]
    iner = ((d^2)/sum((d^2))) * 100
    U = solve(Dr) %*% UDV$u[, 1:r]
    V = solve(Dc) %*% UDV$v[, 1:r]
    D = diagonal(d)
    A = U %*% D
    B = V%*% D
    sf = rowSums(A^2)
    cf = solve(diagonal(sf)) %*% (A^2)
    sc = rowSums(B^2)
    cc = solve(diagonal(sc)) %*% (B^2)
    A = U %*% D^alfa
    B = V %*% D^(1 - alfa)
    afc = list()
    afc$RowCoordinates = A
    afc$ColCoordinates = B
    afc$RowContributions = cf
    afc$ColContributions = cc
    afc$Inertia = iner
    afc$Eigenvalues = d^2
    class(afc) = "afc.sol"
    return(afc)
}
