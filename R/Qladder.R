#' Ladder-of-powers quantile-normal plots
#'
#' Produces a 3x3 matrix of quantile-normal plots for the
#' nine transformations used by Stata's qladder command.
#'
#' @param x Numeric vector. Must contain only positive finite values.
#' @param main Optional overall title.
#' @param digits Number of digits used for numerical labels.
#' @param cex Point size.
#' @param cex.axis Size of axis labels.
#' @param cex.lab Size of axis titles.
#' @param cex.main Size of individual plot titles.
#'
#' @return Invisibly returns a data frame containing the
#' transformed variables.
#'
#' @examples
#' qladder(mtcars$mpg)
#'
#' @export

qladder <- function(
    x,
    main = NULL,
    digits = 2,
    cex = 0.65,
    cex.axis = 0.8,
    cex.lab = 0.8,
    cex.main = 0.9
) {

  # ============================================================
  # 1. VALIDATION
  # ============================================================

  if (!is.numeric(x)) {
    stop("'x' must be numeric.")
  }

  if (length(x) == 0) {
    stop("'x' contains no observations.")
  }

  if (!is.numeric(digits) ||
      length(digits) != 1 ||
      !is.finite(digits) ||
      digits < 0 ||
      digits != floor(digits)) {

    stop("'digits' must be a non-negative integer.")
  }

  # Remove NA, NaN and Inf
  x <- x[is.finite(x)]

  if (length(x) < 3) {
    stop(
      "'x' must contain at least 3 finite observations."
    )
  }

  # ============================================================
  # 2. POSITIVITY REQUIREMENT
  # ============================================================

  if (any(x <= 0)) {
    stop(
      "'x' must contain only positive values. ",
      "The ladder of powers requires x > 0."
    )
  }

  # ============================================================
  # 3. LADDER OF POWERS
  # ============================================================

  transformations <- list(

    "Inverse cubic" =
      1 / x^3,

    "Inverse square" =
      1 / x^2,

    "Inverse" =
      1 / x,

    "Inverse root" =
      1 / sqrt(x),

    "Log" =
      log(x),

    "Square root" =
      sqrt(x),

    "Identity" =
      x,

    "Square" =
      x^2,

    "Cube" =
      x^3
  )

  # ============================================================
  # 4. CHECK TRANSFORMATIONS
  # ============================================================

  for (i in seq_along(transformations)) {

    z <- transformations[[i]]

    if (any(!is.finite(z))) {

      stop(
        "The transformation '",
        names(transformations)[i],
        "' produced non-finite values."
      )
    }
  }

  # ============================================================
  # 5. SAVE GRAPHICAL PARAMETERS
  # ============================================================

  old_par <- par(no.readonly = TRUE)

  on.exit(
    par(old_par),
    add = TRUE
  )

  # ============================================================
  # 6. GRAPH MATRIX
  # ============================================================

  par(
    mfrow = c(3, 3),

    # Bottom, left, top, right
    mar = c(4.5, 4.5, 3.0, 2.0),

    # Outside the complete 3x3 matrix
    oma = c(
      1.5,
      1.5,
      if (is.null(main)) 1 else 3,
      1.5
    )
  )

  # ============================================================
  # 7. Q-Q PLOTS
  # ============================================================

  for (nm in names(transformations)) {

    z <- transformations[[nm]]

    qqnorm(
      z,

      main = nm,

      pch = 16,

      cex = cex,

      cex.axis = cex.axis,

      cex.lab = cex.lab,

      cex.main = cex.main,

      xlab = "Normal quantiles",

      ylab = "Sample quantiles"
    )

    qqline(
      z,
      lwd = 2
    )
  }

  # ============================================================
  # 8. OVERALL TITLE
  # ============================================================

  if (!is.null(main)) {

    mtext(
      main,
      outer = TRUE,
      line = 1,
      font = 2,
      cex = 1.1
    )
  }

  # ============================================================
  # 9. RETURN TRANSFORMED VARIABLES
  # ============================================================

  result <- data.frame(

    invcube =
      transformations[["Inverse cubic"]],

    invsq =
      transformations[["Inverse square"]],

    inv =
      transformations[["Inverse"]],

    invsqrt =
      transformations[["Inverse root"]],

    log =
      transformations[["Log"]],

    sqrt =
      transformations[["Square root"]],

    ident =
      transformations[["Identity"]],

    square =
      transformations[["Square"]],

    cube =
      transformations[["Cube"]]
  )

  invisible(result)
}
