\name{summary.tabu}
\alias{summary.tabu}
\title{
 R Based Tabu Search Summary Function
}
\description{
Summarizes the results of a tabu search optimization run.
}
\usage{
\method{summary}{tabu}(object, verbose = FALSE, ...)
}

\arguments{
  \item{object}{a tabu object.}
  \item{verbose}{if true, the optimal configuration(s) will be printed.}
  \item{\dots}{other options (ignored).}
}
\examples{
# A simple example

evaluateSimple <- function(th)return(1)
result <- tabuSearch(size = 20, iters = 100, objFunc = evaluateSimple) 

summary(result)
summary(result, verbose = TRUE)
}


