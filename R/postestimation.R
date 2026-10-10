# Internal utilities shared by the public post-estimation commands.
.sta_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install the '", package, "' package to use this calculation.", call. = FALSE)
}

.sta_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level))
    stop("level must be one finite number.", call. = FALSE)
  if (level > 1) level <- level / 100
  if (level <= 0 || level >= 1) stop("level must lie strictly between 0 and 1 (or 0 and 100).", call. = FALSE)
  level
}

.sta_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
}

.sta_mode <- function(x) {
  u <- unique(x[!is.na(x)])
  if (!length(u)) stop("Cannot summarize an entirely missing variable.", call. = FALSE)
  u[which.max(tabulate(match(x, u)))]
}

.sta_binary_mean <- function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else .sta_mode(x)

.sta_coef <- function(model) {
  b <- stats::coef(model)
  if (is.matrix(b)) {
    nm <- as.vector(t(outer(rownames(b), colnames(b), paste, sep = ":")))
    b <- stats::setNames(as.vector(t(b)), nm)
  }
  if (inherits(model, "polr")) b <- c(b, model$zeta)
  if (inherits(model, "negbin")) b <- c(b, lnalpha = -log(model$theta))
  if (inherits(model, "fixest") && !is.null(model$theta)) b <- c(b, lnalpha = -log(model$theta))
  b
}

.sta_cov <- function(model, vcov = NULL, df = NULL, vce = NULL, cluster = NULL,
                     hc_type = NULL, correction = NULL) {
  supplied <- inherits(model, "ext_summary")
  if (supplied) {
    summary <- model
    model <- summary$model
    if (is.null(vcov) && is.null(vce) && is.null(cluster) && is.null(hc_type) && is.null(correction)) {
      return(list(model = model, V = summary$vcov,
                  df = if (is.null(df)) summary$inference_df else df, label = summary$vce))
    }
  }
  packages <- c(mlogit="mlogit", clogit="survival", multinom="nnet", polr="MASS",
                negbin="MASS", ivreg="ivreg", fixest="fixest")
  for (class in names(packages)) if (inherits(model,class)) .sta_require(packages[[class]])
  if (inherits(vcov, "ext_summary")) {
    if (!identical(model, vcov$model)) stop("The covariance summary belongs to a different fitted model.", call. = FALSE)
    return(list(model = model, V = vcov$vcov,
                df = if (is.null(df)) vcov$inference_df else df, label = vcov$vce))
  }
  if (!is.null(vcov)) {
    if (!is.null(vce) || !is.null(cluster) || !is.null(hc_type) || !is.null(correction))
      stop("Use vcov OR the vce/cluster/hc_type/correction arguments.", call. = FALSE)
    return(list(model = model, V = vcov, df = if (is.null(df)) Inf else df,
                label = if (is.matrix(vcov)) "Custom" else if (is.function(vcov)) "Custom function" else as.character(vcov)[1]))
  }
  compatible <- inherits(model, c("lm", "glm", "multinom", "polr", "clogit", "mlogit", "ivreg", "fixest"))
  if (compatible) {
    args <- list(model = model, test = "none")
    if (!is.null(vce)) args$vce <- vce
    if (!is.null(cluster)) args$cluster <- cluster
    if (!is.null(hc_type)) args$hc_type <- hc_type
    if (!is.null(correction)) args$correction <- correction
    if (!is.null(df)) args$df <- df
    # Use precisely the same covariance and inference convention as ext_summary.
    result <- tryCatch({
      utils::capture.output(summary <- do.call(ext_summary, args))
      summary
    }, error = identity)
    if (!inherits(result, "error"))
      return(list(model = model, V = result$vcov, df = result$inference_df, label = result$vce))
    if (!is.null(vce) || !is.null(cluster) || !is.null(hc_type) || !is.null(correction)) stop(result)
    # Other links/classes may be supported by marginaleffects, but not ext_summary.
    warning("ext_summary covariance unavailable: ", conditionMessage(result),
            " Using the model's native covariance.", call. = FALSE)
  } else if (!is.null(vce) || !is.null(cluster) || !is.null(hc_type) || !is.null(correction)) {
    stop("For this model, supply vcov explicitly instead of vce options.", call. = FALSE)
  }
  list(model = model, V = stats::vcov(model),
       df = if (is.null(df)) Inf else df, label = "Model-based")
}

.sta_validate_V <- function(V, b) {
  if (!is.matrix(V) || !is.numeric(V) || !identical(dim(V), c(length(b), length(b))) ||
      any(!is.finite(V))) stop("Supply a finite square covariance matrix for all parameters.", call. = FALSE)
  if (!is.null(rownames(V)) && !is.null(names(b))) {
    if (!setequal(rownames(V), names(b)) || !setequal(colnames(V), names(b)))
      stop("Covariance names do not match the parameter names.", call. = FALSE)
    V <- V[names(b), names(b), drop = FALSE]
  }
  if (!isTRUE(all.equal(V, t(V), tolerance = 1e-8, check.attributes = FALSE)))
    stop("The covariance matrix must be symmetric.", call. = FALSE)
  if (min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < -1e-8 * max(1, max(abs(V))))
    stop("The covariance matrix is not positive semidefinite.", call. = FALSE)
  V
}

.sta_delta <- function(fun, b, V, level, df = Inf, null = 0, numderiv = "Richardson") {
  .sta_require("numDeriv")
  if (any(!is.finite(b))) stop("Non-estimable parameters: remove aliased/unidentified coefficients before using the delta method.", call. = FALSE)
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0)
    stop("df must be positive or Inf.", call. = FALSE)
  V <- .sta_validate_V(V, b)
  est <- fun(b)
  if (!is.numeric(est) || !length(est) || any(!is.finite(est)))
    stop("The estimand must return a finite numeric vector.", call. = FALSE)
  if (!is.numeric(null) || !length(null) || !length(null) %in% c(1L, length(est)) || any(!is.finite(null)))
    stop("null must be finite and have length one or the number of estimates.", call. = FALSE)
  J <- numDeriv::jacobian(fun, b, method = numderiv)
  Vout <- J %*% V %*% t(J)
  se <- sqrt(pmax(0, diag(Vout)))
  z <- (as.numeric(est) - null) / se
  p <- if (is.finite(df)) 2 * stats::pt(-abs(z), df) else 2 * stats::pnorm(-abs(z))
  critical <- if (is.finite(df)) stats::qt((1 + level)/2, df) else stats::qnorm((1 + level)/2)
  table <- data.frame(estimate = as.numeric(est), std.error = se, statistic = z, p.value = p,
                      conf.low = est - critical * se, conf.high = est + critical * se,
                      row.names = NULL, check.names = FALSE)
  list(table = table, vcov = Vout, jacobian = J, labels = names(est))
}

.sta_result <- function(table, title, estimate_label, level, df, digits, metadata = list(),
                        notes = character(), covariance = NULL, raw = NULL, class = "statator_margins") {
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits != as.integer(digits) || digits < 0 || digits > 12)
    stop("digits must be an integer from 0 to 12.", call. = FALSE)
  result <- structure(list(table = table, title = title, estimate_label = estimate_label,
                           level = level, df = df, digits = as.integer(digits),
                           metadata = metadata, notes = unique(notes), vcov = covariance, raw = raw),
                      class = c(class, "statator_table"))
  print(result)
  invisible(result)
}

#' Print a StatatoR Post-estimation Table
#' @param x A post-estimation result.
#' @param digits Decimal places for display only.
#' @param ... Reserved for printing methods.
#' @return The result, invisibly.
#' @export
print.statator_table <- function(x, digits = x$digits, ...) {
  fmt <- function(z, d = digits) {
    vapply(z, function(v) {
      if (is.na(v)) return(".")
      if (is.infinite(v)) return(as.character(v))
      if ((v != 0 && abs(v) < 10^(-d)) || abs(v) >= 1e8) formatC(v, format = "e", digits = 2)
      else formatC(v, format = "f", digits = d)
    }, character(1))
  }
  cat("\n", x$title, "\n\n", sep = "")
  for (name in names(x$metadata))
    cat(sprintf("%-24s = %s\n", name, paste(x$metadata[[name]], collapse = ", ")))
  table <- x$table
  numeric <- c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  identifiers <- setdiff(names(table), numeric)
  display <- lapply(table[identifiers], as.character)
  labels <- c(Outcome = "Outcome", Term = "Variable", Comparison = "Comparison",
              Alternative = "Alternative", Case = "Case", Scenario = "Scenario")
  if (length(display)) names(display) <- ifelse(names(display) %in% names(labels), labels[names(display)], names(display))
  display[[x$estimate_label]] <- fmt(table$estimate)
  display[["Std. Error"]] <- fmt(table$std.error)
  display[[if (is.finite(x$df)) "t" else "z"]] <- fmt(table$statistic, 2)
  display[[if (is.finite(x$df)) "P>|t|" else "P>|z|"]] <- fmt(table$p.value, 3)
  display[["Sig."]] <- vapply(table$p.value, function(p) {
    if (is.na(p)) "" else as.character(stats::symnum(p, corr = FALSE, na = FALSE,
      cutpoints = c(0, .001, .01, .05, .1, 1), symbols = c("***", "**", "*", ".", " ")))
  }, character(1))
  display[["CI Lower"]] <- fmt(table$conf.low)
  display[["CI Upper"]] <- fmt(table$conf.high)
  display <- lapply(display, function(v) {v[is.na(v)] <- "."; v})
  widths <- vapply(seq_along(display), function(j) max(nchar(c(names(display)[j], display[[j]])), 5L), integer(1))
  width <- sum(widths) + 3L * (length(widths) - 1L)
  rule <- strrep("-", width)
  cat("\n", rule, "\n", sep = "")
  cat(paste(vapply(seq_along(display), function(j) sprintf("%*s", widths[j], names(display)[j]), character(1)), collapse = "   "), "\n", sep = "")
  cat(rule, "\n", sep = "")
  for (i in seq_len(nrow(table))) {
    if (i > 1L && "Outcome" %in% names(table) && !identical(table$Outcome[i], table$Outcome[i-1L])) cat("\n")
    cat(paste(vapply(seq_along(display), function(j) sprintf("%*s", widths[j], display[[j]][i]), character(1)), collapse = "   "), "\n", sep = "")
  }
  cat(rule, "\n", sep = "")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat(sprintf("Confidence level = %.1f%%\n", 100 * x$level))
  if (length(x$notes)) for (note in x$notes) cat("Note: ", note, "\n", sep = "")
  cat("\n")
  invisible(x)
}

#' Extract a Post-estimation Table
#' @param x A StatatoR post-estimation result.
#' @param row.names Optional row names.
#' @param optional Ignored; retained for the data.frame generic.
#' @param ... Additional arguments passed to as.data.frame.
#' @return A data.frame with unrounded results.
#' @export
as.data.frame.statator_table <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$table, row.names = row.names, optional = optional, ...)
}

#' Nonlinear Combinations of Model Parameters
#'
#' @description General multivariate delta-method estimates, with a formatted table.
#' @param model A fitted model, an ext_summary result, or a named numeric parameter vector.
#' @param expression A function of the named parameter vector returning one or more
#' finite numeric estimates. Names on its return value label the rows.
#' @param vcov Covariance matrix on the original parameter scale. Required for vectors.
#' An ext_summary object may also supply the covariance for its own fitted model.
#' @param null Null values for the transformed estimates, scalar or one per estimate.
#' @param level Confidence level, proportion or percentage.
#' @param df Reference degrees of freedom; Inf gives normal inference.
#' @param digits Decimal places for display only.
#' @param vce,cluster,hc_type,correction Covariance options passed to ext_summary.
#' @param numderiv Numerical Jacobian method passed to numDeriv; default Richardson.
#' @return Invisibly, a statator_table object with table, vcov and Jacobian.
#' @details The transformation must be differentiable near the fitted parameters.
#' Wald intervals are on the transformed scale and may extend beyond natural bounds.
#' This command does not re-estimate the model or check identification of the original fit.
#' @examples
#' m <- glm(am ~ wt, data = mtcars, family = binomial())
#' nlcom(m, function(b) c(OR_wt = exp(b["wt"])), null = 1)
#' @export
nlcom <- function(model, expression, vcov = NULL, null = 0, level = .95,
                  df = NULL, digits = 6, vce = NULL, cluster = NULL,
                  hc_type = NULL, correction = NULL, numderiv = "Richardson") {
  level <- .sta_level(level)
  if (!is.function(expression)) stop("expression must be a function of the parameter vector.", call. = FALSE)
  if (is.numeric(model) && is.null(dim(model))) {
    if (is.null(names(model)) || anyDuplicated(names(model)) || any(!is.finite(model)))
      stop("Supply a finite, uniquely named parameter vector.", call. = FALSE)
    if (is.null(vcov)) stop("A covariance matrix is required for a parameter vector.", call. = FALSE)
    cov <- list(V = vcov, df = if (is.null(df)) Inf else df, label = "Custom")
    b <- model
  } else {
    cov <- .sta_cov(model, vcov, df, vce, cluster, hc_type, correction)
    b <- .sta_coef(cov$model)
  }
  delta <- .sta_delta(expression, b, cov$V, level, cov$df, null, numderiv)
  table <- delta$table
  table$Term <- if (is.null(delta$labels)) paste0("Expression ", seq_len(nrow(table))) else delta$labels
  out <- .sta_result(table, "Nonlinear combinations of parameters", "Estimate", level, cov$df,
                    digits, metadata = list(VCE = cov$label, "Null hypothesis" = null),
                    covariance = delta$vcov, class = "statator_nlcom")
  out$jacobian <- delta$jacobian
  invisible(out)
}
