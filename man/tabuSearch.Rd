\name{tabuSearch}
\alias{tabuSearch}
\title{
R based tabu search algorithm for binary configurations
}
\description{
A tabu search algorithm for optimizing binary strings. It takes a user defined objective function and reports the best binary configuration found throughout the search i.e. the one with the highest objective function value. The algorithm can be used for variable selection. 

The results can be plotted and summarized using \code{plot.tabu} and  \code{summary.tabu}.  }
\usage{
tabuSearch(size = 10, iters = 100, objFunc = NULL, config = NULL, 
neigh = size, listSize = 9, nRestarts = 10, repeatAll = 1, verbose = FALSE)
}

\arguments{
  \item{size}{
the length of the binary configuration.
}
  \item{iters}{
the number of iterations in the preliminary search of the algorithm.
}
  \item{objFunc}{
a user supplied method that evaluates the objective function for a given binary string. The objective function is required to take as an argument a vector of zeros and ones. 
}
  \item{config}{
a starting configuration. 
}
  \item{neigh}{
a number of neighbour configurations to check at each iteration. The default is all, which is the length of the string. If \code{neigh} < \code{size}, the neighbours are chosen at random. 
}
  \item{listSize}{tabu list size. 
}
  \item{nRestarts}{
the maximum number of restarts in the intensification stage of the algorithm.}
  \item{repeatAll}{
the number of times to repeat the search.
}
  \item{verbose}{
if true, the name of the current stage of the algorithm is printed e.g. preliminary stage, intensification stage, diversification stage.
}
}

\value{
A tabu object which is a list giving:

\item{configKeep}{a matrix of configurations chosen at each iteration of the algorithm.}
\item{eUtilityKeep}{value of the objective function for the above configurations.}
\item{iters}{the number of iterations in the preliminary search of the algorithm.}
\item{neigh}{a number of neighbour configurations checked at each iteration.}
\item{listSize}{tabu list size.}
\item{repeatAll}{the number of times the search was repeated.}
}
\references{
Glover, F., 1977. Heuristics for integer programming using surrogate constraints.
Decision Sciences 8, 156-166.



Glover, F., 1986. Future paths for integer programming and links to artificial
intelligence. Computers and Operations Research 13, 533-549.


Fouskakis, D., Draper, D., 2002. Stochastic optimization: a review. International Statistical Review 70, 315-349.
}
\author{K. Domijan}


\seealso{
\code{\link{plot.tabu}} \code{\link{summary.tabu}}
}
\examples{

# A simple example


evaluateSimple <- function(th)return(1)

result <- tabuSearch(size = 20, iters = 100, objFunc = evaluateSimple) 


\dontrun{
# simulate 10-d data: 150 samples from 3 bivariate normals and 8 noise variables. 
# Variable selection should recover the first two variables




library(MASS)
NF <- 10
G <- 3
NTR <- 50
NTE <- 50

muA <- c(1,3)
SigmaA <- matrix(c(0.2, 0.04, 0.04, 0.2), 2, 2)
muB <- c(1.2,1)
SigmaB <- matrix(c(0.1, -0.06, 0.004, 0.2), 2, 2)
muC <- c(3,2)
SigmaC <- matrix(c(0.2, 0.004, 0.004, 0.2), 2, 2)

train <- rbind(mvrnorm(NTR, muA, SigmaA), mvrnorm(NTR, muB, SigmaB), mvrnorm(NTR, muC, SigmaC))
test <- rbind(mvrnorm(NTE, muA, SigmaA), mvrnorm(NTE, muB, SigmaB), mvrnorm(NTE, muC, SigmaC))

train <- cbind(train, matrix(runif(G * NTR * (NF - 2), 0, 4), nrow = G * NTR, ncol = (NF-2)))
test <- cbind(test, matrix(runif(G * NTE * (NF - 2), 0, 4), nrow = G * NTE, ncol = (NF-2)))

wtr <-  as.factor(rep(1:G, each = NTR))
wte <-  as.factor(rep(1:G, each = NTE))
pairs(train, col = wtr)




library(e1071)

evaluate <- function(th){ 
    if (sum(th) == 0)return(0)             
    model <- svm(train[ ,th==1], wtr , gamma = 0.1)
    pred <- predict(model, test[ ,th==1])
    csRate <- sum(pred == wte)/NTE 
    penalty <- (NF - sum(th))/NF 
    return(csRate + penalty)
}  
 
res <- tabuSearch(size = NF, iters = 50, objFunc = evaluate, config = matrix(1,1,NF),
listSize = 5, nRestarts = 4) 


plot(res)
plot(res, "tracePlot")

summary(res, verbose = TRUE)
}
}

