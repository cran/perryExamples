# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Ridge regression with penalty parameter selection
#'
#' Fit ridge regression models and select the penalty parameter by estimating
#' the respective prediction error via (repeated) \eqn{K}-fold
#' cross-validation, (repeated) random splitting (also known as random
#' subsampling or Monte Carlo cross-validation), or the bootstrap.
#'
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param lambda  a numeric vector of non-negative values to be used as penalty
#' parameter.
#' @param standardize  a logical indicating whether the predictor variables
#' should be standardized to have unit variance (the default is \code{TRUE}).
#' @param intercept  a logical indicating whether a constant term should be
#' included in the model (the default is \code{TRUE}).
#' @param splits  an object giving data splits to be used for prediction error
#' estimation (see \code{\link[perry]{perryTuning}}).
#' @param cost  a cost function measuring prediction loss (see
#' \code{\link[perry]{perryTuning}} for some requirements).  The
#' default is to use the root mean squared prediction error (see
#' \code{\link[perry]{cost}}).
#' @param selectBest,seFactor  arguments specifying a criterion for selecting
#' the best model (see \code{\link[perry]{perryTuning}}).  The default is to
#' use a one-standard-error rule.
#' @param ncores,cl  arguments for parallel computing (see
#' \code{\link[perry]{perryTuning}}).
#' @param seed  optional initial seed for the random number generator (see
#' \code{\link{.Random.seed}} and \code{\link[perry]{perryTuning}}).
#' @param \dots  for \code{ridge}, additional arguments to be passed to the
#' prediction loss function \code{cost}.  For \code{ridge.fit}, additional
#' arguments are currently ignored.
#'
#' @return
#' For \code{ridge}, an object of class \code{"perryTuning"}, see
#' \code{\link[perry]{perryTuning}}).  It contains information on the
#' prediction error criterion, and includes the final model with the optimal
#' tuning paramter as component \code{finalModel}.
#'
#' For \code{ridge.fit}, an object of class \code{ridge} with the following
#' components:
#' \describe{
#'   \item{\code{lambda}}{a numeric vector containing the values of the penalty
#'   parameter.}
#'   \item{\code{coefficients}}{a numeric vector or matrix containing the
#'   coefficient estimates.}
#'   \item{\code{fitted.values}}{a numeric vector or matrix containing the
#'   fitted values.}
#'   \item{\code{residuals}}{a numeric vector or matrix containing the
#'   residuals.}
#'   \item{\code{standardize}}{a logical indicating whether the predictor
#'   variables were standardized to have unit variance.}
#'   \item{\code{intercept}}{a logical indicating whether the model includes a
#'   constant term.}
#'   \item{\code{muX}}{a numeric vector containing the means of the predictors.}
#'   \item{\code{sigmaX}}{a numeric vector containing the standard deviations
#'   of the predictors.}
#'   \item{\code{muY}}{numeric; the mean of the response.}
#'   \item{\code{call}}{the matched function call.}
#' }
#'
#' @author Andreas Alfons
#'
#' @references
#' Hoerl, A.E. and Kennard, R.W. (1970) Ridge regression: biased estimation for
#' nonorthogonal problems.  \emph{Technometrics}, \bold{12}(1), 55--67.
#'
#' @seealso
#' \code{\link[perry]{perryTuning}}
#'
#' @example inst/doc/examples/example-ridge.R
#'
#' @keywords regression
#'
#' @export

ridge <- function(x, y, lambda, standardize = TRUE, intercept = TRUE,
                  splits = foldControl(), cost = rmspe,
                  selectBest = c("hastie", "min"), seFactor = 1,
                  ncores = 1, cl = NULL, seed = NULL, ...) {
  # initializations
  if(!is.numeric(lambda) || length(lambda) == 0 || any(!is.finite(lambda))) {
    stop("missing or invalid value of 'lambda'")
  }
  if(any(negative <- lambda < 0)) {
    lambda[negative] <- 0
    warning("negative value for 'lambda', using no penalization")
  }
  lambda <- sort.int(unique(lambda), decreasing=TRUE)
  selectBest <- match.arg(selectBest)
#   call <- call("ridge.fit", lambda=lambda, intercept=intercept,
#                standardize=standardize)
  args <- list(lambda=lambda, intercept=intercept,
               standardize=standardize)
  # estimate prediction error for all values of penalty parameter
#   pe <- perryFit(call, x=x, y=y, splits=splits, cost=cost, costArgs=list(...),
#                  envir=parent.frame(), ncores=ncores, cl=cl, seed=seed)
  pe <- perryFit(ridge.fit, x=x, y=y, args=args, splits=splits, cost=cost,
                 costArgs=list(...), envir=parent.frame(), ncores=ncores,
                 cl=cl, seed=seed)
  # select optimal value of penalty parameter by reshaping results
  pe <- perryReshape(pe, tuning=list(lambda=lambda), selectBest=selectBest,
                     seFactor=seFactor)
  # add final model and return object
  pe$finalModel <- ridge.fit(x, y, lambda=lambda[pe$best], intercept=intercept,
                             standardize=standardize)
  pe
}


#' @rdname ridge
#' @export
#' @import quantreg

ridge.fit <- function(x, y, lambda, standardize = TRUE,
                      intercept = TRUE, ...) {
  # initializations
  matchedCall <- match.call()
  n <- length(y)
  x <- addColnames(as.matrix(x))
  d <- dim(x)
  if(!isTRUE(n == d[1])) stop(sprintf("'x' must have %d rows", n))
  if(!is.numeric(lambda) || length(lambda) == 0 || any(!is.finite(lambda))) {
    stop("missing or invalid value of 'lambda'")
  }
  if(any(negative <- lambda < 0)) {
    lambda[negative] <- 0
    warning("negative value for 'lambda', using no penalization")
  }
  if(length(lambda) > 1) names(lambda) <- seq_along(lambda)
  standardize <- isTRUE(standardize)
  intercept <- isTRUE(intercept)
  # center the data
  if(intercept) {
    muX <- colMeans(x)
    muY <- mean(y)
    xs <- sweep(x, 2, muX, check.margin=FALSE)  # sweep out column centers
    z <- y - muY
  } else {
    muX <- rep.int(0, d[2])
    muY <- 0
    xs <- x
    z <- y
  }
  # scale the predictors
  if(standardize) {
    f <- function(v) sqrt(sum(v^2) / max(1, length(v)-1))
    sigmaX <- apply(xs, 2, f)
    xs <- sweep(xs, 2, sigmaX, "/", check.margin=FALSE)# sweep out column scales
  } else sigmaX <- rep.int(1, d[2])
  # compute ridge solution
  SVD <- svd(xs)
  sv <- SVD$d  # singular values
  a <- sv * drop(crossprod(SVD$u, z)) / (sv^2 + rep(lambda, each=length(sv)))
  dim(a) <- c(length(sv), length(lambda))
  coef <- SVD$v %*% a
  if(length(lambda) == 1) {
    coef <- drop(coef)
    names(coef) <- colnames(x)
  } else dimnames(coef) <- list(colnames(x), names(lambda))
  # back-transform coefficients
  coef <- backtransform(coef, muX=muX, sigmaX=sigmaX, muY=muY,
                        intercept=intercept)
  # compute fitted values and residuals
  fitted <- if(intercept) addIntercept(x) %*% coef else x %*% coef
  fitted <- drop(fitted)
  residuals <- y - fitted
  # construct return object
  fit <- list(lambda=lambda, coefficients=coef, fitted.values=fitted,
              residuals=residuals, standardize=standardize,
              intercept=intercept, muX=muX, sigmaX=sigmaX, muY=muY,
              call=matchedCall)
  class(fit) <- "ridge"
  fit
}
