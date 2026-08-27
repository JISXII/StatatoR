#' Ladder-of-powers quantile-normal plots
#'
#' Produces nine quantile-normal plots corresponding to the
#' ladder of powers used by Stata's qladder command.
#'
#' @param x Numeric vector. Must contain only positive values.
#' @param main Optional overall title.
#' @param digits Number of digits used for numerical labels.
#'
#' @return Invisibly returns a data frame containing the nine
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

  # Box-Cox / ladder transformations require positive values
  if (any(x <= 0)) {
    stop(
      "'x' must contain only positive values. ",
      "The ladder of powers requires x > 0."
    )
  }

  # ============================================================
  # 2. TRANSFORMATIONS
  # ============================================================

  transformations <- list(

    "Inverse cubic" = 1 / x^3,

    "Inverse square" = 1 / x^2,

    "Inverse" = 1 / x,

    "Inverse root" = 1 / sqrt(x),

    "Log" = log(x),

    "Square root" = sqrt(x),

    "Identity" = x,

    "Square" = x^2,

    "Cube" = x^3
  )

  # ============================================================
  # 3. SAVE GRAPHICAL PARAMETERS
  # ============================================================

  old_par <- par(no.readonly = TRUE)

  # Restore graphical parameters when function exits
  on.exit(
    par(old_par),
    add = TRUE
  )

  # ============================================================
  # 4. GRAPHICAL LAYOUT
  # ============================================================

  par(
    mfrow = c(3, 3),
    mar = c(3.5, 3.5, 2.5, 1),
    oma = c(0, 0, ifelse(is.null(main), 0, 2), 0)
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
      col = "blue",
      cex = 0.7,
      xlab = "Normal quantiles",
      ylab = "Sample quantiles"
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
      line = 0.5,
      font = 2,
      cex = 1.2
    )
  }

  # ============================================================
  # 7. RETURN DATA
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
