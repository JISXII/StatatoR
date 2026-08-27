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
#' @param col_pts Color for the points in the Q-Q plot.
#' @param col_line Color for the reference line.
#' @param cex_base Global multiplier for all text sizes.
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
    cex.lab = 1.1,
    cex.main = 0.9,
    col_pts = "darkblue",
    col_line = "red",
    cex_base = 1
) {
  
  # Capture the variable name for the global X-axis
  var_name <- deparse(substitute(x))
  
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
  # 3. LADDER OF POWERS (Stata Exact Order)
  # ============================================================
  
  transformations <- list(
    "cubic" = x^3,
    "square" = x^2,
    "identity" = x,
    "sqrt" = sqrt(x),
    "log" = log(x),
    "1/sqrt" = 1 / sqrt(x),
    "inverse" = 1 / x,
    "1/square" = 1 / x^2,
    "1/cubic" = 1 / x^3
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
    # Reduce inner margins to make it clean like Stata
    mar = c(2.5, 2.5, 2.5, 1.0),
    # Add outer margins for the global labels
    oma = c(4, 4, if (is.null(main)) 3 else 4, 1)
  )
  
  # ============================================================
  # 7. Q-Q PLOTS
  # ============================================================
  
  for (nm in names(transformations)) {
    
    z <- transformations[[nm]]
    
    # Omit inner xlab and ylab to place them globally
    qqnorm(
      z,
      main = nm,
      pch = 16,
      col = col_pts,
      cex = cex * cex_base,
      cex.axis = cex.axis * cex_base,
      cex.main = cex.main * cex_base,
      xlab = "",
      ylab = ""
    )
    
    # 45-degree reference line
    qqline(
      z,
      col = col_line,
      lwd = 2
    )
  }
  
  # ============================================================
  # 8. OVERALL TITLE & GLOBAL AXES
  # ============================================================
  
  if (!is.null(main)) {
    mtext(
      main,
      outer = TRUE,
      line = 1.5,
      font = 2,
      cex = 1.2 * cex_base
    )
  } else {
    mtext(
      "Quantile-Normal plots by transformation",
      outer = TRUE,
      line = 1,
      font = 2,
      cex = 1.2 * cex_base
    )
  }
  
  # Global Y-axis
  mtext(
    "Sample quantiles",
    side = 2,
    outer = TRUE,
    line = 2.2,
    cex = cex.lab * cex_base
  )
  
  # Global X-axis (Variable name)
  mtext(
    var_name,
    side = 1,
    outer = TRUE,
    line = 2.2,
    cex = cex.lab * cex_base
  )
  
  # ============================================================
  # 9. RETURN TRANSFORMED VARIABLES
  # ============================================================
  
  result <- data.frame(
    cube = transformations[["cubic"]],
    square = transformations[["square"]],
    ident = transformations[["identity"]],
    sqrt = transformations[["sqrt"]],
    log = transformations[["log"]],
    invsqrt = transformations[["1/sqrt"]],
    inv = transformations[["inverse"]],
    invsq = transformations[["1/square"]],
    invcube = transformations[["1/cubic"]]
  )
  
  invisible(result)
}
