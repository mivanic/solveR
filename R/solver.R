#' solver: A package solving a system of equations.
#'
#' The solver package provides two important functions: rsolve (which solves a system) and errors (which evaluates error terms)
#'
#' @docType package
#' @name solver
NULL

error <- function(formula, variables) {
    leftSide = eval(formula[[2]], as.list(variables))
    rightSide = eval(formula[[3]], as.list(variables))
    return(leftSide - rightSide)
}

#' Errors of a list of formulas.
#'
#' \code{errors} returns all differences between the left and right sides of a list of formula (errors).
#'
#' @param formulas A list of formulas
#' @param endogenous A list of variables to be used in evaluation (variables may be scalars, vectors or matrices)
#' @param exogenous A list of parameters to be used in evaluation
#' @return A list of real numbers on the same skeleton as the endogenous variables
#' @export
#' @examples
#' errors(formulas =
#'           list(price = price ~ 0.3 * variablecost + coefficient1,
#'                variablecost  = variablecost ~ 0.5 * price + coefficient2),
#'        endogenous = list(price = 1, variablecost=0.4),
#'        exogenous = list(coefficient1=0.56, coefficient2=0.7))

errors <- function(endogenous, exogenous, formulas) {
    e <- sapply(formulas, function(formula) {
        error(formula = formula, variables = c(endogenous, exogenous))
    })
    return(utils::relist(e, endogenous))
}

jacobian <- function(formulas, endogenous, exogenous, pert = 0.001) {

    if (is.null(endogenous))
        stop("There needs to be a vector/list of endogenous variables")

    init = unlist(errors(formulas = formulas, endogenous = endogenous, exogenous = exogenous))
    pert = (rep(stats::runif(1), length(unlist(endogenous))) - 0.5)/100

    sp = unlist(endogenous) %*% t(rep(1, length(unlist(endogenous))))

    diag(sp) = diag(sp) + pert
    rownames(sp) = names(unlist(endogenous))

    updated = apply(sp, MARGIN = 2, function(end) unlist(errors(endogenous = utils::relist(flesh=end, skeleton=endogenous), exogenous = exogenous, formulas = formulas)))

    toRet = (updated - init)/t(pert %*% t(rep(1, length(unlist(endogenous)))))
    return(toRet)
}

doIteration <- function(endogenous, exogenous, formulas) {
    j = jacobian(endogenous = endogenous, exogenous = exogenous, formulas = formulas)
    result = unlist(endogenous) - 0.5 * c(base::solve(j) %*% unlist(errors(formulas = formulas, endogenous = endogenous, exogenous = exogenous)))
    return(utils::relist(flesh = result, skeleton = endogenous))
}

#' Solve a system of equations.
#'
#' \code{solveSystem} returns a set of endogenous variables that solve the system.
#'
#' @param lowerBounds A list of lower bounds for endogenous variables
#' @param upperBounds A list of upper bounds for endogenous variables
#' @param maxIterations A maximum number of times the Newton algorithm should be applied (default value = 100)
#' @param maxError Maximum total absolute error of all equations (default value = 5e-9)
#' @inheritParams errors
#' @return A vector of real numbers
#' @export
#' @examples
#' solveSystem(formulas =
#'           list(price = price ~ 0.3 * variablecost + coefficient1,
#'                variablecost  = variablecost ~ 0.5 * price + coefficient2),
#'        endogenous = list(price = 1, variablecost=0.4),
#'        exogenous = list(coefficient1=0.56, coefficient2=0.7))
#' solveSystem(formulas =
#'           list(price = price ~ 0.3 * variablecost + coefficient1,
#'                variablecost  = variablecost ~ 0.5 * price + coefficient2),
#'        endogenous = list(price = c(wheat=1, rice=4), variablecost=c(wheat=0.4, rice=0.2)),
#'        exogenous = list(coefficient1=c(wheat=0.56, rice=4), coefficient2=c(wehat=0.7, rice=0.9)))

solveSystem <- function(formulas, endogenous, exogenous, lowerBounds = NULL, upperBounds = NULL, maxIterations = 100, maxError = 5e-09) {

    for (it in 1:maxIterations) {
        currentErrors = errors(formulas = formulas, endogenous = endogenous, exogenous = exogenous)
        if (length(unlist(currentErrors)) == 0) {
            currentError = 0
        } else {
            currentError = sum(abs(unlist(currentErrors)))
        }
        if (currentError < maxError)
            break
        endogenous = doIteration(endogenous = endogenous, exogenous = exogenous, formulas = formulas)


        unlistedEndogenous = unlist(endogenous)
        unlistedLowerBounds = unlist(lowerBounds)
        unlistedUpperBounds = unlist(upperBounds)

        # if the endogenous variable falls below the lower bound, move it above
        unlistedEndogenous[unlistedEndogenous < unlistedLowerBounds[names(unlistedEndogenous)]] = stats::runif(length(unlistedEndogenous[unlistedEndogenous <
            unlistedLowerBounds[names(unlistedEndogenous)]]))/10 + unlistedLowerBounds[names(unlistedEndogenous[unlistedEndogenous < unlistedLowerBounds[names(unlistedEndogenous)]])]

        # if the endogenous variable exceeds the upperbound, move it above
        unlistedEndogenous[unlistedEndogenous > unlistedUpperBounds[names(unlistedEndogenous)]] = -stats::runif(length(unlistedEndogenous[unlistedEndogenous >
            unlistedUpperBounds[names(unlistedEndogenous)]]))/10 + unlistedUpperBounds[names(unlistedEndogenous[unlistedEndogenous > unlistedUpperBounds[names(unlistedEndogenous)]])]

        endogenous = utils::relist(flesh = unlistedEndogenous, skeleton = endogenous)

    }
    if (currentError >= maxError) {
        stop(paste("Error too large (over", maxError, ") after", maxIterations, "iterations"))
    }

    return(list(results = endogenous, totalError = currentError, iterations = it))
}

#' Presolve a system--i.e., find which variables may be calculated directly and return a set of solved variables and the remaining equations.
#'
#' \code{presolveSystem} returns all differences between the left and right sides of a list of formula (errors).
#'
#' @param formulas A list of formulas
#' @return A list of remaining formulas and solved variables
#' @export
#' @examples
#' errors(formulas =
#'           list(price = price ~ 0.3 * variablecost + coefficient1,
#'                variablecost  = variablecost ~ 0.5 * price + coefficient2),
#'        endogenous = list(price = 1, variablecost=0.4),
#'        exogenous = list(coefficient1=0.56, coefficient2=0.7))

