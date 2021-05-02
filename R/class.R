#' Survival Summary Object.
#'
#' @slot Marginal Per arm statistics.
#' @slot Contrasts Contrasts of the marginal statistics.
#' @name SurvSummary-class
#' @rdname SurvSummary-class
#' @exportClass SurvSummary

setClass(
  Class = "SurvSummary",
  representation = representation(
   Marginal = "data.frame",
   Contrasts = "data.frame"
  )
)

# -----------------------------------------------------------------------------
# Print Method
# -----------------------------------------------------------------------------

#' Print Method for SurvSummary Object.
#'
#' Print method for objects of class \code{SurvSummary}.
#'
#' @param x An object of class \code{SurvSummary}.
#' @param ... Unused.
#' @export

print.SurvSummary <- function (x, ...) {
  
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

# -----------------------------------------------------------------------------
# Show Method
# -----------------------------------------------------------------------------

#' Show Method for SurvSummary Object
#'
#' @param object An object of class \code{SurvSummary}.
#' @rdname fit-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "SurvSummary"),
  definition = function (object) {print.SurvSummary(x = object)}
)

