# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' LAD-lasso with penalty parameter selection
#'
#' Fit LAD-lasso models and select the penalty parameter by estimating the
#' respective prediction error via (repeated) \eqn{K}-fold cross-validation,
#' (repeated) random splitting (also known as random subsampling or Monte Carlo
#' cross-validation), or the bootstrap.
#'
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param lambda  for \code{ladlasso}, a numeric vector of non-negative values
#' to be used as penalty parameter.  For \code{ladlasso.fit}, a single
#' non-negative value to be used as penalty parameter.
#' @param standardize  a logical indicating whether the predictor variables
#' should be standardized to have unit MAD (the default is \code{TRUE}).
#' @param intercept  a logical indicating whether a constant term should be
#' included in the model (the default is \code{TRUE}).
#' @param splits  an object giving data splits to be used for prediction error
#' estimation (see \code{\link[perry]{perryTuning}}).
#' @param cost  a cost function measuring prediction loss (see
#' \code{\link[perry]{perryTuning}} for some requirements).  The
#' default is to use the mean absolute prediction error (see
#' \code{\link[perry]{cost}}).
#' @param selectBest,seFactor  arguments specifying a criterion for selecting
#' the best model (see \code{\link[perry]{perryTuning}}).  The default is to
#' use a one-standard-error rule.
#' @param ncores,cl  arguments for parallel computing (see
#' \code{\link[perry]{perryTuning}}).
#' @param seed  optional initial seed for the random number generator (see
#' \code{\link{.Random.seed}} and \code{\link[perry]{perryTuning}}).
#' @param \dots  for \code{ladlasso}, additional arguments to be passed to the
#' prediction loss function \code{cost}.  For \code{ladlasso.fit}, additional
#' arguments to be passed to \code{\link[quantreg]{rq.fit.lasso}}.
#'
#' @return
#' For \code{ladlasso}, an object of class \code{"perryTuning"}, see
#' \code{\link[perry]{perryTuning}}).  It contains information on the
#' prediction error criterion, and includes the final model with the optimal
#' tuning paramter as component \code{finalModel}.
#'
#' For \code{ladlasso.fit}, an object of class \code{ladlasso} with the
#' following components:
#' \describe{
#'   \item{\code{lambda}}{numeric; the value of the penalty parameter.}
#'   \item{\code{coefficients}}{a numeric vector containing the coefficient
#'   estimates.}
#'   \item{\code{fitted.values}}{a numeric vector containing the fitted values.}
#'   \item{\code{residuals}}{a numeric vector containing the residuals.}
#'   \item{\code{standardize}}{a logical indicating whether the predictor
#'   variables were standardized to have unit MAD.}
#'   \item{\code{intercept}}{a logical indicating whether the model includes a
#'   constant term.}
#'   \item{\code{muX}}{a numeric vector containing the medians of the
#'   predictors.}
#'   \item{\code{sigmaX}}{a numeric vector containing the MADs of the
#'   predictors.}
#'   \item{\code{muY}}{numeric; the median of the response.}
#'   \item{\code{call}}{the matched function call.}
#' }
#'
#' @author Andreas Alfons
#'
#' @references
#' Wang, H., Li, G. and Jiang, G. (2007) Robust regression shrinkage and
#' consistent variable selection through the LAD-lasso.  \emph{Journal of
#' Business & Economic Statistics}, \bold{25}(3), 347--355.
#'
#' @seealso
#' \code{\link[perry]{perryTuning}}, \code{\link[quantreg]{rq.fit.lasso}}
#'
#' @example inst/doc/examples/example-ladlasso.R
#'
#' @keywords regression robust
#'
#' @export

ladlasso <- function(x, y, lambda, standardize = TRUE, intercept = TRUE,
                     splits = foldControl(), cost = mape,
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
#   call <- call("ladlasso.fit", intercept=intercept, standardize=standardize)
  args <- list(intercept=intercept, standardize=standardize)
  # estimate prediction error for each value of penalty parameter
  # and add final model
#   perryTuning(call, x=x, y=y, tuning=list(lambda=lambda), splits=splits,
#               cost=cost, costArgs=list(...), selectBest=selectBest,
#               seFactor=seFactor, final=TRUE, envir=parent.frame(),
#               ncores=ncores, cl=cl, seed=seed)
  perryTuning(ladlasso.fit, x=x, y=y, args=args, tuning=list(lambda=lambda),
              splits=splits, cost=cost, costArgs=list(...),
              selectBest=selectBest, seFactor=seFactor, final=TRUE,
              envir=parent.frame(), ncores=ncores, cl=cl, seed=seed)
}


#' @rdname ladlasso
#' @export
#' @import quantreg

ladlasso.fit <- function(x, y, lambda, standardize = TRUE,
                         intercept = TRUE, ...) {
  # initializations
  matchedCall <- match.call()
  n <- length(y)
  x <- addColnames(as.matrix(x))
  d <- dim(x)
  if(!isTRUE(n == d[1])) stop(sprintf("'x' must have %d rows", n))
  lambda <- rep_len(lambda, length.out=1)
  if(!is.numeric(lambda) || !is.finite(lambda)) {
    stop("missing or invalid value of 'lambda'")
  }
  if(lambda < 0) {
    lambda <- 0
    warning("negative value for 'lambda', using no penalization")
  }
  standardize <- isTRUE(standardize)
  intercept <- isTRUE(intercept)
  # center the data
  if(intercept) {
    muX <- apply(x, 2, median)
    muY <- median(y)
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
    sigmaX <- apply(xs, 2, mad, center=0)
    # fallback mode in case of zero MAD
    tooSmall <- which(sigmaX <= .Machine$double.eps)
    if(length(tooSmall) > 0) {
      # center with mean if requested
      if(intercept) {
        # sweep out column means
        muX[tooSmall] <- colMeans(x[, tooSmall, drop=FALSE])
        xs[, tooSmall] <- sweep(x[, tooSmall, drop=FALSE], 2, muX[tooSmall],
                                check.margin=FALSE)
      }
      # standardize with standard deviation
      f <- function(v) sqrt(sum(v^2) / max(1, length(v)-1))
      sigmaX[tooSmall] <- apply(xs[, tooSmall, drop=FALSE], 2, f)
    }
    # sweep out column scales
    xs <- sweep(xs, 2, sigmaX, "/", check.margin=FALSE)
  } else sigmaX <- rep.int(1, d[2])
  # compute LAD-lasso solution
  rqlasso <- function(..., tau) rq.fit.lasso(..., tau=0.5)
  fit <- rqlasso(xs, z, lambda=rep_len(lambda, d[2]), ...)
  # back-transform coefficients
  coef <- backtransform(coef(fit), muX=muX, sigmaX=sigmaX, muY=muY,
                        intercept=intercept)
  # compute fitted values and residuals
  residuals <- drop(residuals(fit))
  fitted <- y - residuals
  # construct return object
  fit <- list(lambda=lambda, coefficients=coef, fitted.values=fitted,
              residuals=residuals, standardize=standardize,
              intercept=intercept, muX=muX, sigmaX=sigmaX, muY=muY,
              call=matchedCall)
  class(fit) <- "ladlasso"
  fit
}
