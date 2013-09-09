\name{Env}
\alias{Env}
\title{
  Ecological Factors in Farm Management
}
\description{
  The farmsEnv data frame has 20 rows and 4 columns. The rows are farms
  on the Dutch island of Terschelling and the columns are factors
  describing the management of grassland.
}
\usage{
data(Env)
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Mois}{five levels of soil moisture}
    \item{Manag}{Grassland management type (SF = standard,
      BF = biological, HF = hobby farming,
      NM = nature conservation).
    }
    \item{Use}{Grassland use (U1 = it exists production, U2 =
      intermediate,U3 = grazing).
    }
    \item{Manure}{
      Manure usage (C0,C1,C2,C3 and C4)
    }
  }
}
\source{
  J.C. Gower and D.J. Hand (1996) \emph{Biplots}. Chapman & Hall, Table 4.6.

  Quoted as from:\cr
  R.H.G. Jongman, C.J.F. ter Braak and O.F.R. van Tongeren (1987)
  \emph{Data Analysis in Community and Landscape Ecology.}
  PUDOC, Wageningen.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\examples{
data(Env)
}
\keyword{datasets}
