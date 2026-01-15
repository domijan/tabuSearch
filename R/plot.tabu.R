#' R Based Tabu Search Plot Function
#'
#' @param x a tabu object.
#' @param type "tracePlot" or default.
#' @param ... options directly passed to the plot function.
#'
#' @returns NULL
#' @export
#'@description
#' Plots features of an optimization run of the tabu search algorithm for binary strings. The default plots show: (a) the number of times each element of the string was set to one over the search, (b) frequency of moves for each element of the string over the serach, (c) the number of ones in the chosen configuration at each iteration, (d) the objective function value of the current configuration at each iteration of the algorithm. The "tracePlot" shows the current configurations for all interations.
#'
#' @examples
#' # A simple example
#'
#' evaluateSimple <- function(th)return(1)
#' result <- tabuSearch(size = 20, iters = 100, objFunc = evaluateSimple)
#'
#' plot(result)
#' plot(result, "tracePlot")

plot.tabu <-
function (x, type = "default", ...)
{
    tabuObject <- x

    if (type == "tracePlot"){
        image(x=seq(1,dim(tabuObject$configKeep)[2], by=1), y = seq(1,dim(tabuObject$configKeep)[1], by=1), z=t(-tabuObject$configKeep),
        xlab = "variable", ylab = "iterations", col = gray(0:8 / 8),   useRaster  = TRUE)
        }

    else if (type == "default"){
        def.par <- par(no.readonly = TRUE)
        nf<-layout(matrix(c(1,2,3,4), 2, 2))
        nSelect <- colSums(tabuObject$configKeep)
        plot(nSelect, xlab = "variable", ylab ="", main = "No of times selected", type = "l")
        mostFrequentMoves <- apply(tabuObject$configKeep, 2, function(y)sum(diff(y) != 0))
        plot(mostFrequentMoves, xlab = "variable", ylab ="", main = "Most frequent moves", type = "l")
        nVars <- rowSums(tabuObject$configKeep)
        plot(nVars, xlab = "iterations", ylab ="", main = "Sum of included variables", type = "l")
        plot(tabuObject$eUtilityKeep, xlab = "iterations", ylab ="", main = "Objective Function", type = "l")
        par(def.par)#- reset to default
}
    else  stop("error: Plot type not supported for a tabu object")
    invisible()
}
