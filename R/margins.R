#' Predictions, Marginal Effects and Discrete Changes
#'
#' @description A compact, formatted interface to marginaleffects for ordinary
#' regressions, with separate conditional/choice-model calculations. No predictor,
#' outcome, alternative or number of categories is hard-coded.
#' @param model A fitted model or an ext_summary result. The latter preserves its VCE.
#' @param dydx Predictor names for derivatives. NULL requests predictions. "all"
#' selects predictors; do not combine with change.
#' @param change Named list of two values per predictor, e.g. list(treated=c(0,1)),
#' for discrete changes (second value minus first). Mutually exclusive with dydx.
#' @param at Named list of predictor values defining scenarios. Multiple vectors
#' generate their Cartesian product. At observed values, other variables retain
#' their observed distribution; with atmeans=TRUE they are summarized instead.
#' @param atmeans Summarize numeric predictors at their means and categorical
#' predictors at their modes. Numeric 0/1 predictors also use their means.
#' @param newdata Optional prediction data. Choice-model data must retain case and
#' alternative indices, or supply a dfidx object. For complex formulas supply the
#' original raw variables, not only transformed model-matrix columns.
#' @param average Average over prediction observations/cases; default TRUE.
#' @param by Optional grouping variable names. For mlogit, groups must be constant
#' within each choice case. Combined with scenario and outcome grouping.
#' @param outcome Optional outcome names to retain in the returned table.
#' @param type Prediction scale. Defaults to probs for multinom/polr, pc1 for
#' clogit, probabilities for mlogit, and response otherwise. Clogit additionally
#' supports pu0: inverse-logit of the index with fixed effect zero, not pc1.
#' @param vcov,vce,cluster,hc_type,correction Covariance inputs. Default uses
#' ext_summary conventions where supported; otherwise native model covariance.
#' vcov can also be FALSE (ordinary backend only), a matrix, or a backend request.
#' @param level Confidence level, proportion or percentage.
#' @param df Reference degrees of freedom; defaults to the covariance convention.
#' @param digits Decimal places for printing; calculations are unrounded.
#' @param numderiv Delta-method differentiation: richardson (default) or a method
#' supported by marginaleffects. Custom choice calculations use Richardson.
#' @param eps Positive numeric step for continuous predictor differentiation.
#' NULL uses a scale-dependent step. Check sensitivity for ill-scaled variables.
#' @param slope dydx, eyex, eydx or dyex, delegated to marginaleffects.
#' Choice-model calculations currently support dydx only.
#' @param variable_type Optional named character vector, case or alternative,
#' overriding automatic within-case predictor classification for choice models.
#' @param ... Additional marginaleffects arguments, such as wts, hypothesis or
#' comparison. Reserved wrapper arguments cannot be duplicated. Not used by the
#' custom clogit/mlogit engine, which rejects them rather than ignoring them.
#' @return Invisibly, a statator_margins object with table, metadata, notes and raw
#' backend output. Extract unrounded results with result$table or as.data.frame().
#' @details Ordinary models use marginaleffects, including its treatment of
#' binary/categorical predictors as discrete changes when supplied to dydx.
#' Use change to state the two comparison values explicitly. The choice engine
#' treats numeric dydx predictors as continuous, including numeric 0/1 variables.
#' For mlogit, changing an alternative attribute produces own and cross effects;
#' changing a case attribute changes all alternatives in that case together.
#' Only deterministic ordinary/nested mlogit predictions are supported by this
#' engine; random-coefficient, probit and heteroscedastic choice models are rejected.
#' Clogit pc1 supports predictions only, requires one chosen alternative per stratum,
#' and cannot be atmeans. Use pu0 for clogit derivatives or discrete changes.
#' Custom choice calculations require unweighted models, a matrix covariance, and
#' no offsets and unconstrained mlogit fits; they do not silently approximate unsupported specifications.
#' @examples
#' m <- glm(am ~ wt, data=mtcars, family=binomial())
#' margins(m)
#' margins(m, dydx="wt", atmeans=TRUE)
#' margins(m, at=list(wt=c(2,3,4)))
#' @export
margins <- function(model, dydx = NULL, change = NULL, at = NULL, atmeans = FALSE,
                    newdata = NULL, average = TRUE, by = NULL, outcome = NULL,
                    type = NULL, vcov = NULL, vce = NULL, cluster = NULL,
                    hc_type = NULL, correction = NULL, level = .95, df = NULL,
                    digits = 6, numderiv = "richardson", eps = NULL,
                    slope = "dydx", variable_type = NULL, ...) {
  .sta_require("marginaleffects")
  level <- .sta_level(level)
  .sta_flag(atmeans, "atmeans"); .sta_flag(average, "average")
  if (!is.null(dydx) && !is.null(change)) stop("Use dydx OR change, not both.", call. = FALSE)
  if (!is.null(dydx) && (!is.character(dydx) || !length(dydx) || anyNA(dydx) || anyDuplicated(dydx)))
    stop("dydx must contain distinct predictor names, or 'all'.", call. = FALSE)
  if (!is.null(change) && (!is.list(change) || !length(change) || is.null(names(change)) ||
      any(!nzchar(names(change))) || anyDuplicated(names(change)) ||
      any(lengths(change) != 2L) || any(vapply(change, anyNA, logical(1)))))
    stop("change must be a named list of two nonmissing values per predictor.", call. = FALSE)
  if (!is.null(at) && (!is.list(at) || !length(at) || is.null(names(at)) ||
      any(!nzchar(names(at))) || anyDuplicated(names(at)) || any(lengths(at) == 0L) ||
      any(vapply(at, anyNA, logical(1)))))
    stop("at must be a named, nonempty list of nonmissing scenario values.", call. = FALSE)
  if (!is.null(by) && (!is.character(by) || anyNA(by) || anyDuplicated(by)))
    stop("by must contain distinct grouping variable names.", call. = FALSE)
  if (!is.null(eps) && (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0))
    stop("eps must be one positive finite number.", call. = FALSE)
  cov <- .sta_cov(model, vcov, df, vce, cluster, hc_type, correction)
  model <- cov$model
  quantity <- if (!is.null(change)) "Discrete changes" else if (!is.null(dydx)) "Marginal effects" else "Predictions"
  if (is.null(type)) type <- if (inherits(model, c("multinom", "polr"))) "probs" else
    if (inherits(model, "clogit")) "pc1" else if (inherits(model, "mlogit")) "probabilities" else "response"
  if (inherits(model, c("clogit", "mlogit"))) {
    if (!identical(tolower(numderiv), "richardson")) stop("Custom choice margins support numderiv='richardson' only.", call. = FALSE)
    if (length(list(...))) stop("Additional backend arguments are not supported for clogit/mlogit margins.", call. = FALSE)
    if (slope != "dydx") stop("Choice-model margins currently support slope='dydx' only.", call. = FALSE)
    result <- .sta_choice_margins(model, cov, dydx, change, at, atmeans, newdata,
                                  average, by, type, level, eps, variable_type)
    table <- result$table; raw <- NULL; notes <- result$notes; outcov <- result$vcov
  } else {
    if (!is.null(variable_type)) stop("variable_type applies to choice models only.", call. = FALSE)
    data <- newdata
    if (atmeans || !is.null(at)) {
      gridargs <- c(list(FUN_binary = .sta_binary_mean), at)
      if (is.null(data)) gridargs$model <- model else gridargs$newdata <- data
      if (!atmeans) gridargs$grid_type <- "counterfactual"
      data <- do.call(marginaleffects::datagrid, gridargs)
    }
    groups <- unique(c(by, names(at)))
    if (inherits(model, c("multinom", "polr")) && length(groups)) groups <- unique(c("group", groups))
    args <- list(model = model, newdata = data, type = type, vcov = cov$V,
                 conf_level = level, df = cov$df, numderiv = numderiv)
    # Count-model predictions do not directly depend on dispersion. Preserve the
    # joint beta covariance while removing the auxiliary-only row for the backend.
    if (is.matrix(args$vcov) && "lnalpha" %in% rownames(args$vcov) &&
        !"lnalpha" %in% names(stats::coef(model)))
      args$vcov <- args$vcov[rownames(args$vcov) != "lnalpha", colnames(args$vcov) != "lnalpha", drop = FALSE]
    if (!average) args$by <- FALSE
    else if (length(groups)) args$by <- groups
    if (!is.null(change)) {
      args$variables <- change
      fun <- if (average) marginaleffects::avg_comparisons else marginaleffects::comparisons
    } else if (!is.null(dydx)) {
      args$variables <- if (identical(dydx, "all")) NULL else dydx
      args$slope <- slope
      if (!is.null(eps)) args$eps <- eps
      fun <- if (average) marginaleffects::avg_slopes else marginaleffects::slopes
    } else fun <- if (average) marginaleffects::avg_predictions else marginaleffects::predictions
    extra <- list(...)
    if (length(extra) && (is.null(names(extra)) || any(!nzchar(names(extra))) || anyDuplicated(names(extra))))
      stop("Additional arguments must be uniquely named.", call. = FALSE)
    if (any(names(extra) %in% names(args))) stop("A reserved backend argument was supplied twice.", call. = FALSE)
    raw <- do.call(fun, c(args, extra))
    all <- as.data.frame(raw)
    numeric <- c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
    for (name in setdiff(numeric, names(all))) all[[name]] <- NA_real_
    ids <- intersect(unique(c("group", "term", "contrast", groups, if (!average) "rowid")), names(all))
    table <- all[, c(ids, numeric), drop = FALSE]
    names(table)[match(c("group", "term", "contrast", "rowid"), names(table), nomatch = 0L)] <-
      c(group = "Outcome", term = "Term", contrast = "Comparison", rowid = "Case")[c("group", "term", "contrast", "rowid") %in% names(table)]
    notes <- if (!is.null(dydx)) "Binary/categorical predictors follow marginaleffects discrete-change conventions; use change for explicit comparisons." else character()
    outcov <- NULL
  }
  if (!is.null(outcome)) {
    if (!"Outcome" %in% names(table) || !all(outcome %in% as.character(table$Outcome)))
      stop("outcome must match outcome names in the result.", call. = FALSE)
    keep <- as.character(table$Outcome) %in% outcome
    table <- table[keep, , drop = FALSE]
    if (!is.null(outcov)) outcov <- outcov[keep, keep, drop = FALSE]
  }
  if ("Outcome" %in% names(table)) {
    ord <- order(table$Outcome)
    table <- table[ord, , drop = FALSE]
    if (!is.null(outcov)) outcov <- outcov[ord, ord, drop = FALSE]
  }
  rownames(table) <- NULL
  probability <- type %in% c("probs", "probabilities", "pc1", "pu0") ||
    (inherits(model,"mlogit") && type=="response") ||
    (inherits(model, "glm") && stats::family(model)$family %in% c("binomial", "quasibinomial") && type == "response") ||
    (inherits(model, "fixest") && ((is.list(model$family) && model$family$family == "binomial") ||
      identical(model$family,"logit")) && type == "response")
  estimate_label <- if (!is.null(change)) "Change" else if (!is.null(dydx))
    switch(slope, dydx = "dy/dx", eyex = "ey/ex", eydx = "ey/dx", dyex = "dy/ex", "Estimate")
    else if (probability) "Probability" else "Prediction"
  .sta_result(table, paste(if (average) "Average" else "Conditional", tolower(quantity)),
              estimate_label, level, cov$df, digits,
              metadata = list(Model = class(model)[1], Type = type, VCE = cov$label,
                              Evaluation = if (atmeans) "At means / categorical modes" else "Observed or supplied values",
                              Scenarios = if (is.null(at)) "None" else paste(names(at), vapply(at, paste, character(1), collapse = "/"), sep = "=")),
              notes = notes, covariance = outcov, raw = raw)
}

#' Plot a Saved Margins Result
#' @param x A result from statator::margins.
#' @param xvar Table column for the horizontal axis. Auto-selects a scenario
#' column, otherwise Outcome, Variable or Case. Specify it when ambiguous.
#' @param group Optional table column identifying separate series.
#' @param facet Optional table column for separate panels.
#' @param ci Include confidence intervals.
#' @param connect Connect ordered numeric x values within a series.
#' @param horizontal Flip coordinate axes.
#' @param title Optional plot title.
#' @param ... Reserved; currently unused.
#' @return A ggplot object, editable using ordinary ggplot2 layers.
#' @examples
#' m <- glm(am ~ wt, mtcars, family=binomial())
#' r <- margins(m, at=list(wt=c(2,3,4)))
#' marginsplot(r, xvar="wt")
#' @export
marginsplot <- function(x, xvar = NULL, group = NULL, facet = NULL,
                        ci = TRUE, connect = TRUE, horizontal = FALSE, title = NULL, ...) {
  .sta_require("ggplot2")
  if (!inherits(x, "statator_margins")) stop("x must be a result from statator::margins().", call. = FALSE)
  .sta_flag(ci, "ci"); .sta_flag(connect, "connect"); .sta_flag(horizontal, "horizontal")
  if (length(list(...))) stop("Add ggplot2 layers to customize the returned plot.", call. = FALSE)
  table <- x$table
  ids <- setdiff(names(table), c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
  if (is.null(xvar)) {
    candidates <- setdiff(ids, c("Outcome", "Term", "Comparison", "Alternative"))
    if (length(candidates) > 1L) stop("Several scenario columns are available; specify xvar.", call. = FALSE)
    xvar <- if (length(candidates)) candidates[1] else if ("Outcome" %in% ids) "Outcome" else if ("Term" %in% ids) "Term" else NULL
  }
  if (is.null(xvar)) {table$Result <- seq_len(nrow(table)); xvar <- "Result"}
  requested <- c(xvar, group, facet)
  if (any(!requested %in% names(table))) stop("Plot variables must name columns in result$table.", call. = FALSE)
  if (any(lengths(list(xvar, group, facet)) > 1L)) stop("Specify one column per plot argument.", call. = FALSE)
  other <- setdiff(ids, requested)
  varying <- other[vapply(table[other], function(v) length(unique(v)) > 1L, logical(1))]
  if (length(varying) && is.null(group)) {
    table$.series <- do.call(interaction, c(table[varying], list(drop = TRUE, sep = " / ")))
    group <- ".series"
  }
  table$.x <- table[[xvar]]
  table$.group <- if (is.null(group)) factor("Estimate") else factor(table[[group]])
  plot <- ggplot2::ggplot(table, ggplot2::aes(x = .data[[".x"]], y = .data[["estimate"]],
                                           colour = .data[[".group"]], group = .data[[".group"]]))
  position <- ggplot2::position_dodge(width = if (is.numeric(table$.x)) 0 else .35)
  if (connect && is.numeric(table$.x)) plot <- plot + ggplot2::geom_line()
  if (ci) plot <- plot + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[["conf.low"]], ymax = .data[["conf.high"]]),
                                              width = 0, position = position, na.rm = TRUE)
  plot <- plot + ggplot2::geom_point(position = position) +
    ggplot2::labs(x = xvar, y = x$estimate_label, colour = if (identical(group, ".series")) "Series" else group,
                  title = if (is.null(title)) x$title else title) + ggplot2::theme_classic()
  if (is.null(group)) plot <- plot + ggplot2::guides(colour = "none")
  if (!is.null(facet)) plot <- plot + ggplot2::facet_wrap(stats::reformulate(facet))
  if (horizontal) plot <- plot + ggplot2::coord_flip()
  plot
}

#' @importFrom rlang .data
NULL
