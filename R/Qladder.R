#' Ladder-of-powers quantile-normal plots
#'
#' Produces the nine quantile-normal plots used by Stata's qladder
#' command.
#'
#' The transformations are:
#' inverse-cubic, inverse-square, inverse, inverse-root,
#' log, square-root, identity, square, and cube.
#'
#' @param x Numeric vector.
#' @param main Optional overall title.
#' @param digits Number of digits used in plot labels.
#'
#' @return Invisibly returns a data frame containing all nine
#' transformed variables.
#'
#' @examples
#' qladder(mtcars$mpg)
#'
#' @export

qladder <- function(x,
                    main = NULL,
                    digits = 2) {

  # ============================================================
  # 1. VALIDATION
  # ============================================================

  if (!is.numeric(x)) {
    stop("'x' must be numeric.")
  }

  if (!is.numeric(digits) ||
      length(digits) != 1 ||
      digits < 0 ||
      digits != floor(digits)) {
    stop("'digits' must be a non-negative integer.")
  }

  # Remove missing and non-finite observations
  x <- x[is.finite(x)]

  if (length(x) < 3) {
    stop("'x' must contain at least 3 finite observations.")
  }

  # All transformations used by qladder require x > 0
  if (any(x <= 0)) {
    stop(
      "'x' must contain only positive values. ",
      "The ladder of powers requires x > 0."
    )
  }

  # ============================================================
  # 2. LADDER OF POWERS
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
  # 3. SAVE CURRENT GRAPHICS SETTINGS
  # ============================================================

  old_par <- par(no.readonly = TRUE)

  on.exit(
    par(old_par),
    add = TRUE
  )

  # ============================================================
  # 4. GRAPH MATRIX
  # ============================================================

  par(
    mfrow = c(3, 3),
    oma = c(0, 0, 3, 0)
  )

  # ============================================================
  # 5. Q-Q PLOTS
  # ============================================================

  for (nm in names(transformations)) {

    z <- transformations[[nm]]

    qqnorm(
      z,
      main = nm,
      pch = 16,
      col = "blue"
    )

    qqline(
      z,
      col = "red",
      lwd = 2
    )
  }

  # ============================================================
  # 6. OVERALL TITLE
  # ============================================================

  if (!is.null(main)) {

    mtext(
      main,
      outer = TRUE,
      line = 1,
      font = 2,
      cex = 1.2
    )
  }

  # ============================================================
  # 7. RETURN TRANSFORMED DATA
  # ============================================================

  result <- data.frame(
    invcube = transformations[["Inverse cubic"]],
    invsq = transformations[["Inverse square"]],
    inv = transformations[["Inverse"]],
    invsqrt = transformations[["Inverse root"]],
    log = transformations[["Log"]],
    sqrt = transformations[["Square root"]],
    ident = transformations[["Identity"]],
    square = transformations[["Square"]],
    cube = transformations[["Cube"]]
  )

  invisible(result)
}
