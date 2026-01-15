#' R Based Tabu Search Summary Function
#'
#'
#' @param object a tabu object.
#' @param verbose if true, the optimal configuration(s) will be printed.
#' @param ... other options (ignored).
#'
#' @returns a matrix of configurations from iterations with maximum value of objective function (i.e. Optimum configuration).
#' @export
#' @description Summarizes the results of a tabu search optimization run.
#' @examples
#' # A simple example
#'
#' evaluateSimple <- function(th)return(1)
#' result <- tabuSearch(size = 20, iters = 100, objFunc = evaluateSimple)
#'
#' summary(result)
#' summary(result, verbose = TRUE)
#'
summary.tabu <-
function (object, verbose = FALSE, ...)
{
    tabuObject <- object

    nVars <- rowSums(tabuObject$configKeep)
    nSelect <- colSums(tabuObject$configKeep)
    uniqueConfig <- dim(unique(tabuObject$configKeep))[1]


    output <- paste("Tabu Settings", "\n",
    "  Type                                       = ", tabuObject$type, "\n",
    "  No of algorithm repeats                    = ", tabuObject$repeatAll, "\n",
    "  No of iterations at each prelim search     = ", tabuObject$iters, "\n",
    "  Total no of iterations                     = ", length(nVars), "\n",
    "  No of unique best configurations           = ", uniqueConfig, "\n",
    "  Tabu list size                             = ", tabuObject$listSize, "\n",
    "  Configuration length                       = ", length(nSelect), "\n",
    "  No of neighbours visited at each iteration = ", tabuObject$neigh, "\n",   sep = "")

    maxObjFunction <- max(tabuObject$eUtilityKeep)
    optimumNVars <- nVars[which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))]
    optimumIteration <- which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))
    optimumConfig <- tabuObject$configKeep[optimumIteration, ]

    optionPart <- paste("Results:", "\n",
    "  Highest value of objective fn    = ", round(maxObjFunction, 5), "\n",
    "  Occurs # of times                = ", length(optimumNVars), "\n",
    "  Optimum number of variables      = ", deparse(optimumNVars),  "\n",  sep = "")

    cat(output)
    cat(optionPart)
    if (verbose) {
        cat(paste("Optimum configuration:", "\n"))
        print(optimumConfig)
        }
}
