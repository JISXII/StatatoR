.sta_choice_data <- function(model, newdata) {
  if (inherits(model, "mlogit")) {
    .sta_require("dfidx")
    if (!is.null(model$rpar) || isTRUE(model$probit) || isTRUE(model$heterosc) ||
        !is.null(model$call$rpar) || isTRUE(model$call$probit) || isTRUE(model$call$heterosc) ||
        isTRUE(model$call$ranked) || isTRUE(model$call$panel))
      stop("Choice margins support deterministic ordinary/nested mlogit, not random/probit/heteroscedastic models.", call. = FALSE)
    if (!inherits(model$model, "dfidx")) stop("Fit mlogit with a retained dfidx model frame.", call. = FALSE)
    if (any(attr(model$coefficients,"fixed") %in% TRUE))
      stop("Choice margins currently require an unconstrained mlogit fit.", call. = FALSE)
    ix <- dfidx::idx(model$model)
    if (ncol(ix) != 2L) stop("Choice margins require exactly case and alternative indices.", call. = FALSE)
    id <- names(ix)[1]; alt <- names(ix)[2]
    response <- all.vars(stats::formula(model))[1]
    data <- if (is.null(newdata)) model$model else newdata
    if (inherits(data, "dfidx")) {
      idx <- dfidx::idx(data)
      data <- as.data.frame(data)
      data$idx <- NULL
      data[[id]] <- idx[[1]]; data[[alt]] <- idx[[2]]
    } else data <- as.data.frame(data)
    if (!all(c(id, alt) %in% names(data)))
      stop("newdata must contain index columns '", id, "' and '", alt, "', or be dfidx.", call. = FALSE)
    if (anyNA(data[, c(id, alt)]) || anyDuplicated(data[, c(id, alt)]))
      stop("Each case/alternative pair must be unique and nonmissing.", call. = FALSE)
    alternatives <- levels(ix[[2]])
    if (is.null(alternatives)) alternatives <- unique(as.character(ix[[2]]))
    if (any(!as.character(data[[alt]]) %in% alternatives)) stop("Unknown alternatives in newdata.", call. = FALSE)
    # predict.mlogit requires complete alternative sets in each scenario.
    sets <- split(as.character(data[[alt]]), data[[id]])
    if (any(!vapply(sets, function(v) setequal(v, alternatives), logical(1))))
      stop("Prediction data must contain all model alternatives in each case. Supply a balanced scenario.", call. = FALSE)
    if (!response %in% names(data)) data[[response]] <- as.character(data[[alt]]) == alternatives[1]
    return(list(data = data, id = id, alt = alt, response = response, alternatives = alternatives))
  }
  .sta_require("insight")
  data <- if (is.null(newdata)) insight::get_data(model) else as.data.frame(newdata)
  if (is.null(data) || !nrow(data)) stop("Cannot recover prediction data; supply newdata explicitly.", call. = FALSE)
  terms <- attr(stats::terms(model), "term.labels")
  strata <- terms[grepl("^(survival::)?strata\\(", terms)]
  if (!length(strata)) stop("Cannot identify clogit strata.", call. = FALSE)
  list(data = data, strata = strata, id = NULL, alt = NULL,
       response = NULL, alternatives = NULL)
}

.sta_mlogit_call <- function(model) {
  # Native update.mlogit evaluates the original call in its caller. Freeze the
  # estimator and normalization so a saved unqualified call is self-contained.
  model$call[[1]] <- quote(mlogit::mlogit)
  base <- .ext_mlogit(model)$base
  if (length(base)) model$call$reflevel <- base[1]
  if (!is.null(model$nests)) model$call$nests <- model$nests
  remove <- intersect(names(model$call), c("choice","shape","varying","sep","alt.var",
    "chid.var","alt.levels","opposite","drop.index","id","ranked"))
  model$call[remove] <- NULL
  model
}

.sta_choice_predict <- function(model, data, info, type) {
  if (inherits(model, "mlogit")) {
    if (!type %in% c("probabilities", "probs", "response"))
      stop("mlogit margins currently return probabilities only.", call. = FALSE)
    data[[info$alt]] <- factor(data[[info$alt]], levels = info$alternatives)
    nd <- dfidx::dfidx(data, shape = "long", choice = info$response, idx = c(info$id, info$alt))
    model <- .sta_mlogit_call(model)
    prob <- stats::predict(model, newdata = nd)
    if (is.null(dim(prob))) prob <- matrix(prob, nrow = 1L, dimnames = list(NULL, info$alternatives))
    cases <- unique(as.character(dfidx::idx(nd)[[1]]))
    if (!is.null(rownames(prob)) && setequal(rownames(prob), cases)) cases <- rownames(prob)
    if (nrow(prob) != length(cases) || is.null(colnames(prob)))
      stop("Unexpected prediction dimensions from mlogit.", call. = FALSE)
    ids <- expand.grid(Case = cases, Outcome = colnames(prob), stringsAsFactors = FALSE)
    for (name in setdiff(names(data), c(info$alt, info$response))) {
      values <- data[[name]][match(ids$Case, as.character(data[[info$id]]))]
      if (!name %in% names(ids)) ids[[name]] <- values
    }
    return(list(estimate = as.numeric(prob), ids = ids))
  }
  if (!type %in% c("pu0", "pc1")) stop("clogit type must be 'pc1' or 'pu0'.", call. = FALSE)
  lp <- as.numeric(stats::predict(model, newdata = data, type = "lp", reference = "zero"))
  if (type == "pu0") prob <- stats::plogis(lp) else {
    env <- new.env(parent = environment(stats::formula(model)))
    env$strata <- survival::strata
    strata <- lapply(info$strata, function(term) eval(parse(text = term)[[1]], data, env))
    groups <- do.call(interaction, c(strata, list(drop = TRUE)))
    exp_lp <- exp(lp - ave(lp, groups, FUN = max))
    prob <- exp_lp / ave(exp_lp, groups, FUN = sum)
  }
  ids <- data
  ids$Case <- seq_len(nrow(data))
  list(estimate = prob, ids = ids)
}

.sta_choice_aggregate <- function(pred, average, by, alternatives = NULL) {
  ids <- pred$ids
  groups <- unique(c(intersect("Outcome", names(ids)), by))
  if (any(!groups %in% names(ids))) stop("Unknown grouping variable.", call. = FALSE)
  if (!average) {
    keys <- ids[, unique(c("Case", groups)), drop = FALSE]
    return(list(estimate = pred$estimate, keys = keys))
  }
  if (!length(groups)) return(list(estimate = mean(pred$estimate), keys = data.frame(.row = 1L)))
  keys <- unique(ids[, groups, drop = FALSE])
  key <- function(d) do.call(paste, c(lapply(d, as.character), sep = "\r"))
  g <- match(key(ids[, groups, drop = FALSE]), key(keys))
  est <- vapply(seq_len(nrow(keys)), function(i) mean(pred$estimate[g == i]), numeric(1))
  list(estimate = est, keys = keys)
}

.sta_mlogit_design <- function(model, data, info) {
  model <- .sta_mlogit_call(model)
  data[[info$alt]] <- factor(data[[info$alt]], levels=info$alternatives)
  nd <- dfidx::dfidx(data, shape="long", choice=info$response, idx=c(info$id, info$alt))
  fit <- stats::update(model, data=nd, start=stats::coef(model, fixed=TRUE), iterlim=0, print.level=0)
  list(X=stats::model.matrix(fit), index=dfidx::idx(fit$model))
}

.sta_choice_margins <- function(model, cov, dydx, change, at, atmeans, newdata,
                                average, by, type, level, eps, variable_type) {
  .sta_require("numDeriv")
  if (!is.matrix(cov$V)) stop("Custom choice margins require a covariance matrix (or an ext_summary result).", call. = FALSE)
  if (!is.null(stats::model.weights(model$model)) || !is.null(stats::model.offset(model$model)))
    stop("Custom choice margins currently require unweighted models without offsets.", call. = FALSE)
  info <- .sta_choice_data(model, newdata)
  data <- info$data
  predictors <- unique(all.vars(stats::formula(model)[[3]]))
  predictors <- intersect(setdiff(predictors, c(info$id, info$alt)), names(data))
  if (inherits(model, "clogit")) {
    # The expanded survival formula includes strata variables; they are not slopes.
    stratum_vars <- unique(unlist(lapply(info$strata, function(s) all.vars(parse(text = s)))))
    predictors <- setdiff(predictors, stratum_vars)
    if (type == "pc1") {
      if (atmeans) stop("At-means pc1 destroys the choice set; use supplied choice-set newdata or type='pu0'.", call. = FALSE)
      if (!is.null(dydx) || !is.null(change))
        stop("clogit pc1 supports predictions only. Use pu0 for its margins, or mlogit for conditional own/cross effects.", call. = FALSE)
      Y <- stats::model.response(model$model)
      grp <- model$model[, grep("strata\\(", names(model$model)), drop = FALSE]
      group <- do.call(interaction, c(grp, list(drop = TRUE)))
      if (any(ave(Y[, ncol(Y)], group, FUN = sum) != 1))
        stop("pc1 is implemented for one chosen alternative per stratum only.", call. = FALSE)
    }
  }
  requested <- if (!is.null(change)) names(change) else if (identical(dydx, "all")) predictors else dydx
  if (any(!c(requested, names(at), by) %in% names(data))) stop("Unknown predictor, scenario or grouping variable.", call. = FALSE)
  if (any(!requested %in% predictors)) stop("Effects must refer to raw predictors, not indices or outcomes.", call. = FALSE)
  if (!is.null(dydx) && any(!vapply(data[requested], is.numeric, logical(1))))
    stop("Use change for nonnumeric predictors.", call. = FALSE)
  types <- stats::setNames(rep("case", length(predictors)), predictors)
  if (inherits(model, "mlogit")) {
    for (name in predictors) {
      groups <- split(data[[name]], data[[info$id]])
      types[name] <- if (all(vapply(groups, function(v) length(unique(v)) <= 1L, logical(1)))) "case" else "alternative"
    }
    if (!is.null(variable_type)) {
      if (is.null(names(variable_type)) || anyDuplicated(names(variable_type)) ||
          any(!names(variable_type) %in% predictors) || any(!variable_type %in% c("case", "alternative")))
        stop("variable_type must be a named vector of 'case' or 'alternative' for predictors.", call. = FALSE)
      types[names(variable_type)] <- variable_type
    }
    for (name in by) {
      if (any(vapply(split(data[[name]], data[[info$id]]), function(v) length(unique(v)) > 1L, logical(1))))
        stop("Choice-model by variables must be constant within each case.", call. = FALSE)
    }
  } else if (!is.null(variable_type)) stop("variable_type applies to mlogit only.", call. = FALSE)
  if (atmeans) {
    if (length(by)) stop("For grouped at-means choice calculations supply explicit newdata scenarios.", call. = FALSE)
    summarize <- function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else .sta_mode(x)
    if (inherits(model, "mlogit")) {
      base <- data[match(info$alternatives, as.character(data[[info$alt]])), , drop = FALSE]
      for (name in predictors) {
        if (types[name] == "case") {
          value <- summarize(data[[name]][!duplicated(data[[info$id]])])
          base[[name]] <- rep(value, nrow(base))
        } else {
          values <- lapply(info$alternatives, function(a) summarize(data[[name]][as.character(data[[info$alt]]) == a]))
          if (is.factor(data[[name]])) base[[name]] <- factor(vapply(values, as.character, character(1)),
            levels=levels(data[[name]]), ordered=is.ordered(data[[name]]))
          else base[[name]] <- unlist(values, use.names=FALSE)
        }
      }
      base[[info$id]] <- 1L
      data <- base
    } else {
      base <- data[1, , drop = FALSE]
      for (name in predictors) base[[name]] <- summarize(data[[name]])
      data <- base
    }
  }
  scenarios <- if (is.null(at)) data.frame(.scenario = 1L) else expand.grid(at, stringsAsFactors = FALSE)
  b <- .sta_coef(model)
  keys <- NULL
  design_cache <- new.env(parent=emptyenv())
  evaluate <- function(b) {
    m <- model
    if (inherits(m, "mlogit")) m$coefficients <- b else m$coefficients <- b
    estimates <- list(); labels <- list(); n <- 0L
    for (s in seq_len(nrow(scenarios))) {
      d <- data
      for (name in names(at)) {
        value <- scenarios[[name]][s]
        if (is.factor(d[[name]])) {
          if (!as.character(value) %in% levels(d[[name]])) stop("Unknown factor level in at.", call. = FALSE)
          d[[name]][] <- as.character(value)
        } else d[[name]][] <- value
      }
      append <- function(pred, term = NULL, alternative = NULL) {
        a <- .sta_choice_aggregate(pred, average, by)
        key <- a$keys
        key$.row <- NULL
        if (!ncol(key)) key <- data.frame(.row = seq_along(a$estimate))
        if (!is.null(term)) key$Term <- term
        if (!is.null(alternative)) key$Alternative <- alternative
        for (name in names(at)) key[[name]] <- scenarios[[name]][s]
        n <<- n + 1L; estimates[[n]] <<- a$estimate; labels[[n]] <<- key
      }
      if (!length(requested)) append(.sta_choice_predict(m, d, info, type)) else {
        for (variable in requested) {
          alternatives <- if (inherits(m, "mlogit") && types[variable] == "alternative") info$alternatives else NA_character_
          for (alternative in alternatives) {
            rows <- if (is.na(alternative)) rep(TRUE, nrow(d)) else as.character(d[[info$alt]]) == alternative
            lo <- hi <- d
            if (!is.null(change)) {
              values <- change[[variable]]
              if (is.factor(d[[variable]]) && any(!as.character(values) %in% levels(d[[variable]]))) stop("Unknown factor level in change.", call. = FALSE)
              lo[[variable]][rows] <- values[1]; hi[[variable]][rows] <- values[2]
              denominator <- 1
            } else {
              step <- if (is.null(eps)) 1e-4 * max(1, stats::sd(info$data[[variable]]), na.rm = TRUE) else eps
              lo[[variable]][rows] <- lo[[variable]][rows] - step
              hi[[variable]][rows] <- hi[[variable]][rows] + step
              denominator <- 2 * step
            }
            cache_key <- paste(s, match(variable, requested), if (is.na(alternative)) 0L else match(alternative, info$alternatives), sep=".")
            if (inherits(m, "clogit") && type == "pu0" && is.null(change)) {
              # Differentiate the design, then use the exact logistic derivative.
              # This avoids a second finite difference of nearly equal probabilities.
              if (!exists(cache_key, design_cache, inherits=FALSE)) {
                Xlo <- stats::model.matrix(m, data = lo)
                Xhi <- stats::model.matrix(m, data = hi)
                assign(cache_key, (Xhi-Xlo)[,names(b),drop=FALSE]/denominator, design_cache)
              }
              derivative_lp <- as.numeric(get(cache_key, design_cache) %*% b)
              high <- .sta_choice_predict(m, d, info, type)
              high$estimate <- high$estimate * (1-high$estimate) * derivative_lp
            } else if (inherits(m,"mlogit") && is.null(m$nests) && is.null(change)) {
              if (!exists(cache_key, design_cache, inherits=FALSE)) {
                Xlo <- .sta_mlogit_design(model,lo,info)
                Xhi <- .sta_mlogit_design(model,hi,info)
                assign(cache_key,list(D=(Xhi$X-Xlo$X)[,names(b),drop=FALSE]/denominator,
                  index=Xhi$index),design_cache)
              }
              cached <- get(cache_key,design_cache)
              derivative_u <- as.numeric(cached$D %*% b)
              high <- .sta_choice_predict(m,d,info,type)
              row <- match(paste(high$ids$Case,high$ids$Outcome,sep="\r"),
                paste(cached$index[[1]],cached$index[[2]],sep="\r"))
              derivative_u <- derivative_u[row]
              P <- high$estimate
              high$estimate <- P*(derivative_u-ave(P*derivative_u,high$ids$Case,FUN=sum))
            } else {
              low <- .sta_choice_predict(m, lo, info, type)
              high <- .sta_choice_predict(m, hi, info, type)
              high$estimate <- (high$estimate-low$estimate)/denominator
            }
            append(high, variable, if (is.na(alternative)) NULL else alternative)
          }
        }
      }
    }
    keys <<- do.call(rbind, labels)
    unlist(estimates, use.names = FALSE)
  }
  delta <- .sta_delta(evaluate, b, cov$V, level, cov$df)
  table <- cbind(keys, delta$table)
  table$.row <- NULL
  notes <- if (inherits(model, "clogit") && type == "pu0")
    "pu0 = inverse-logit(index) with fixed effect zero; this is not a within-choice-set probability." else character()
  if (length(requested) && inherits(model, "mlogit")) notes <- c(notes,
    "Alternative identifies the attribute changed; Outcome identifies the probability affected. Both own and cross effects are reported.",
    paste("Predictor types:", paste(names(types)[names(types) %in% requested], types[names(types) %in% requested], sep = "=", collapse = ", ")))
  if (length(requested) && is.null(change)) notes <- c(notes, "Numeric predictors are differentiated as continuous; check eps sensitivity for unusual scales.")
  list(table = table, vcov = delta$vcov, notes = notes)
}
