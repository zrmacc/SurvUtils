# -----------------------------------------------------------------------------
# One Sample
# -----------------------------------------------------------------------------

#' One Sample Survival Summary Object.
#'
#' @slot CumHaz Cumulative hazard curve.
#' @slot CumHazVar Variance of cumulative hazard.
#' @slot CumHazLower Lower bound of cumulative hazard.
#' @slot CumHazUpper Upper bound of cumulative hazard.
#' @slot NAR Number at risk curve.
#' @slot Surv Survival curve. 
#' @slot SurvVar Variance of survival curve.
#' @slot SurvLower Lower bound of survival curve.
#' @slot SurvUpper Upper bound of survival curve.
#' @slot tmax Max observation time.
#' @name OneSample-class
#' @rdname OneSample-class
#' @exportClass OneSample

setClass(
  Class = "OneSample",
  representation = representation(
    CumHaz = "function",
    CumHazVar = "function",
    CumHazLower = "function",
    CumHazUpper = "function",
    NAR = "function",
    Surv = "function",
    SurvVar = "function",
    SurvLower = "function",
    SurvUpper = "function",
    tmax = "numeric"
  )
)


#' Print Method for OneSample Object.
#'
#' Print method for objects of class \code{OneSample}.
#'
#' @param x An object of class \code{OneSample}.
#' @param ... Unused.
#' @export

print.OneSample <- function (x, ...) {
  
  # Number at-risk.
  nar <- x@NAR(0)
  cat("Number initially at risk: ")
  cat(nar)
  cat("\n")
  
  # Maximum observation time.
  tmax <- x@tmax
  cat("Maximum observation time (tmax): ")
  cat(signif(tmax, digits = 3))
  cat("\n")
  
  # Cumulative hazard at tmax.
  cat("Cumulative hazard at tmax: ")
  cat(signif(x@CumHaz(tmax), digits = 3))
  cat("\n")
  
  # Survival at tmax.
  cat("Survival at tmax: ")
  cat(signif(x@Surv(tmax), digits = 3))
  cat("\n\n")
  
}


#' Show Method for OneSample Object
#'
#' @param object An object of class \code{OneSample}.
#' @rdname OneSample-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "OneSample"),
  definition = function (object) {print.OneSample(x = object)}
)

# -----------------------------------------------------------------------------
# Two Sample.
# -----------------------------------------------------------------------------

#' Two Sample Survival Summary Object.
#'
#' @slot Marginal Per arm statistics.
#' @slot Contrasts Contrasts of the marginal statistics.
#' @name TwoSample-class
#' @rdname TwoSample-class
#' @exportClass TwoSample

setClass(
  Class = "TwoSample",
  representation = representation(
   Marginal = "data.frame",
   Contrasts = "data.frame"
  )
)


#' Print Method for TwoSample Object.
#'
#' Print method for objects of class \code{TwoSample}.
#'
#' @param x An object of class \code{TwoSamples}.
#' @param ... Unused.
#' @export

print.TwoSample <- function (x, ...) {
  
  disp <- function(y) {
    if (is.numeric(y)) {
      out <- signif(y, digits = 3)
    } else {
      out <- y
    }
    return(out)
  }
  
  # Marginal.
  marg <- x@Marginal
  marg[, ] <- lapply(marg, disp)
  cat('Marginal Statistics:\n')
  show(marg)
  cat('\n\n')
  
  # Contrasts.
  contrast <- x@Contrasts
  contrast[, ] <- lapply(contrast, disp)
  cat('Contrasts:\n')
  show(contrast)
  cat('\n\n')

}


#' Show Method for TwoSample Object
#'
#' @param object An object of class \code{TwoSample}.
#' @rdname TwoSample-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "TwoSample"),
  definition = function (object) {print.TwoSample(x = object)}
)

