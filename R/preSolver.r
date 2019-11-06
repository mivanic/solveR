#' Pre solve a system of formulas.
#'
#' \code{presolveSystem} returns a set of remaining formulas and implied constants and functions
#' @param constants (optional) A list of coefficients
#' @return A list of formulas, constants and functions
#' @export
#' @examples
#' presolveSystem(formulas =
#'           list(price = 1 ~ price + coefficient1,
#'                variablecost  = variablecost ~ 0.5 * price + coefficient2), constants=list(coefficient1=3))

presolveSystem <- function(formulas, constants = list()) {
  formulasAfterConstants = evaluateConstants(formulas = formulas, constants =
                                               constants)
  formulasAfterImpliedConstants = do.call(evaluateImpliedConstants, formulasAfterConstants)
  return (do.call(evaluateFunctions, formulasAfterImpliedConstants))
}

#presolveSystem(formulas=list(a=a~b+c, b=a-c~b), constants=list(c=4))

# An internal function to turn a formula into a function
formulaToFunction <- function(formula, constants = list()) {
  toReturn = list()
  allVariables = all.vars(formula[[3]])
  toReturn[[formula[[2]]]] = list(
    `function` = eval(parse(
      text = paste(
        'function(',
        paste(setdiff(allVariables, names(constants)), collapse = ','),
        ')',
        deparse(formula[[3]])
      )
    )),
    representation = paste(formula[[2]], '(',
                           paste(
                             setdiff(allVariables, names(constants)), collapse = ','
                           ),
                           ')')
  )
  return(toReturn)
}

#formulaToFunction(yy~3*log(pp)+oo, constants=list(oo=4))

# An internal function that takes a list of formulas and for marked formulas, it turns them into functions and makes replacements
reduceFormulas = function(formulaList,
                          environmentList,
                          reduceVector) {
  formulaReductionList = Reduce(function(a, e)
    c(a, e), Map(function(y, n) {
      if (length(which(n == reduceVector)) > 0)
        return(formulaToFunction(formula = y, constants = environmentList))
      else{
        return(NULL)
      }
    }, formulaList, names(formulaList)))


  toReturn = list(
    equationList = Filter(function(x)
      ! is.null(x), Map(
        function(y, n) {
          if (length(which(n == reduceVector)) > 0)
            return(NULL)
          else
            return(y)
        }, formulaList, names(formulaList)
      )),
    environmentList = Reduce(
      function(a, e) {
        if (is.null(e))
          return(a)
        else{
          a[[names(e)[1]]] = e[[1]]
          return(a)
        }
      },
      Map(function(y, n) {
        if (length(which(n == reduceVector)) > 0)
          return(formulaToFunction(y)$`function`)
        else{
          return(NULL)
        }
      }, formulaList, names(formulaList)),
      environmentList
    ),
    formulaReductionList = formulaReductionList
  )

  toReturn$equationListTransformed = lapply(toReturn$equationList, function(x) {
    do.call(substitute, list(
      expr = x,
      env = lapply(formulaReductionList, function(y) {
        return(do.call(quote, list(as.formula(
          paste('.', y$representation, sep = '~')
        )[[3]])))
      })
    ))
  })

  return(toReturn)
}

# An internal function that removes formulas that define constants

evaluateConstants = function(formulas,
                             constants = list(),
                             functions = list()) {
  evaluatedOne = F
  for (f in names(formulas)) {
    if (length(setdiff(all.vars(formulas[[f]][[3]]), names(constants))) == 0) {
      # no variable on the right
      constants[[formulas[[f]][[2]]]] = eval(formulas[[f]][[3]], envir =
                                               constants)
      formulas[[f]] = NULL
      evaluatedOne = T
    }
  }
  if (evaluatedOne == T) {
    return(Recall(formulas, constants, functions))
  } else{
    return(list(
      formulas = formulas,
      constants = constants,
      functions = functions
    ))
  }
}

#evaluateConstants(formulas=list(a=w~3, b=u~2+p))

# An internal function that removes formulas that imply  constants (e.g., 6 ~ p - 5 implies that p = 11)

evaluateImpliedConstants = function(formulas,
                                    constants = list(),
                                    functions = list()) {
  evaluatedOne = F
  for (f in names(formulas)) {
    if (length(setdiff(all.vars(formulas[[f]]), names(constants))) == 1) {
      # there is exactly one unknown in formula f
      unknown = setdiff(all.vars(formulas[[f]]), names(constants))
      trySolve = solveSingleFormula(formula = formulas[[f]],
                                    unknown = unknown,
                                    constants = constants)
      if (trySolve$success == T) {
        constants[[unknown]] = trySolve$value
        formulas[[f]] = NULL
        evaluatedOne = T
      }
    }
  }
  if (evaluatedOne == T) {
    return(Recall(formulas, constants, functions))
  } else{
    return(list(
      formulas = formulas,
      constants = constants,
      functions = functions
    ))
  }
}

#evaluateImpliedConstants(formulas=list(a=1+w~3, b=u~2+p))

solveSingleFormula = function(formula, unknown, constants) {
  value1 = 1
  constants[[unknown]] = value1
  error1 = eval(formula[[2]], envir = constants) - eval(formula[[3]], envir = constants)

  value2 = 2
  constants[[unknown]] = value2
  error2 = eval(formula[[2]], envir = constants) - eval(formula[[3]], envir = constants)

  de = (error2 - error1) / (value2 - value1)

  value3 = value2 - error2 / de
  constants[[unknown]] = value3
  error3 = eval(formula[[2]], envir = constants) - eval(formula[[3]], envir = constants)
  if (abs(error3) < 1e-6) {
    return(list(success = T, value = value3))
  } else{
    return(list(success = F))
  }
}

#solveSingleFormula(formula=4~+u+3*y, 'y', list(u=1))

# An internal function that removes formulas that can be turned into functions

evaluateFunctions = function(formulas,
                             constants = list(),
                             functions = list()) {
  evaluatedOne = F
  for (f in names(formulas)) {
    cf = formulas[[f]]
    if (length(cf[[2]]) == 1 &&
        length(all.vars(cf[[2]])) == 1 &&
        is.null(constants[[all.vars(cf[[2]])[[1]]]]) &&
        length(intersect(all.vars(cf[[2]]), all.vars(cf[[3]]))) == 0) {
      reducedFormulas = reduceFormulas(
        formulaList = formulas,
        environmentList = constants,
        reduceVector = c(f)
      )
      evaluatedOne = T
      break
    }
  }
  if (evaluatedOne == T) {
    return(
      Recall(
        formulas = reducedFormulas$equationListTransformed,
        constants = constants,
        functions = c(
          functions,
          Map(
            function(ff)
              ff$`function`,
            reducedFormulas$formulaReductionList
          )
        )
      )
    )
  } else {
    return(list(
      formulas = formulas,
      constants = constants,
      functions = functions
    ))
  }
}
