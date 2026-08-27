#' Generalized Ladder of Powers Histograms (Stata-style)
#'
#' @description
#' Displays a 3x3 grid of histograms for the 9 Tukey ladder of powers transformations,
#' overlaid with a normal density curve. Allows customizing bar colors, line colors,
#' and scaling all font sizes proportionally for presentations or reports.
#'
#' @param x A numeric vector to be transformed and plotted.
#' @param bins Integer. Suggested number of bins for the histograms. Default is 15.
#' @param col_bars Character. Color for the histogram bars. Default is "#b3cde3".
#' @param col_line Character. Color for the overlaid normal curve. Default is "#e41a1c".
#' @param cex_base Numeric. Multiplier for all text sizes (titles, axes, labels). Default is 1.
#'
#' @return A 3x3 base R plot is displayed. No object is returned.
#' @export
gladder <- function(x, bins = 15, col_bars = "#b3cde3", col_line = "#e41a1c", cex_base = 1) {
  
  if (!is.numeric(x)) stop("'x' must be a numeric vector.")
  var_name <- deparse(substitute(x))
  x <- na.omit(x)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(4, 4, 3, 1))
  
  trans_names <- c("cubic", "square", "identity", "square root", "log", 
                   "1/(square root)", "inverse", "1/square", "1/cubic")
  trans_forms <- c("v^3", "v^2", "v", "sqrt(v)", "log(v)", 
                   "1/sqrt(v)", "1/v", "1/(v^2)", "1/(v^3)")
  
  apply_trans <- function(v, type) {
    suppressWarnings(
      switch(type,
             "cubic" = v^3, "square" = v^2, "identity" = v,
             "square root" = sqrt(v), "log" = log(v), "1/(square root)" = 1/sqrt(v),
             "inverse" = 1/v, "1/square" = 1/v^2, "1/cubic" = 1/v^3)
    )
  }
  
  for (i in 1:9) {
    xt <- apply_trans(x, trans_names[i])
    xt <- na.omit(xt)
    xt <- xt[is.finite(xt) & !is.na(xt)]
    
    if (length(xt) < 2 || var(xt) == 0) {
      plot.new()
      title(main = trans_names[i], sub = trans_forms[i], col.main = "darkgray", 
            cex.main = 1 * cex_base, cex.sub = 0.9 * cex_base)
      text(0.5, 0.5, "Invalid Data\n(NAs produced)", col = "red", cex = 0.8 * cex_base)
      next
    }
    
    # Dibujar histograma escalando el título interno y los números de los ejes
    hist(xt, breaks = bins, prob = TRUE, col = col_bars, border = "white",
         main = trans_names[i], xlab = "", ylab = "", 
         cex.main = 1 * cex_base, cex.axis = 0.85 * cex_base)
    
    x_seq <- seq(min(xt), max(xt), length.out = 100)
    y_seq <- dnorm(x_seq, mean = mean(xt), sd = sd(xt))
    lines(x_seq, y_seq, col = col_line, lwd = 2)
  }
  
  # Añadir etiquetas globales por fuera de la grilla 3x3 escaladas por cex_base
  mtext("Histograms by transformation", outer = TRUE, line = 1, cex = 1.2 * cex_base, font = 2)
  mtext("Density", side = 2, outer = TRUE, line = 2.2, cex = 1.1 * cex_base)
  mtext(var_name, side = 1, outer = TRUE, line = 2.2, cex = 1.1 * cex_base)
}
