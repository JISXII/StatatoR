.sta_quadratic <- function(d, V, tolerance = 1e-9) {
  V <- (V + t(V))/2
  eig <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(V)) || !length(d) ||
      min(abs(eig)) <= tolerance * max(abs(eig)))
    return(list(statistic = NA_real_, p = NA_real_, status = "Singular covariance difference"))
  q <- as.numeric(crossprod(d, solve(V, d)))
  invalid <- q < 0 || min(eig) < -tolerance * max(abs(eig))
  list(statistic = q, p = if (invalid) NA_real_ else stats::pchisq(q, length(d), lower.tail = FALSE),
       status = if (q < 0) "Negative chi2: not interpretable" else if (invalid)
         "Indefinite covariance difference: not interpretable" else "OK")
}

.sta_iia_row <- function(method, omitted, d, V) {
  q <- .sta_quadratic(d, V)
  data.frame(Test = method, Omitted = omitted, statistic = q$statistic,
              df = length(d), p.value = q$p, Status = q$status, stringsAsFactors = FALSE)
}

.sta_mnl_vector <- function(m) {
  b <- stats::coef(m)
  if (is.matrix(b)) return(stats::setNames(as.vector(t(b)),
    as.vector(t(outer(rownames(b), colnames(b), paste, sep = ":")))))
  stats::setNames(b, paste(m$lev[2], names(b), sep = ":"))
}

.sta_mnl_vcov <- function(m) {
  V <- stats::vcov(m)
  dimnames(V) <- rep(list(names(.sta_mnl_vector(m))), 2)
  V
}

.sta_mnl_scores <- function(m) {
  X <- stats::model.matrix(m$terms, m$model, contrasts.arg = m$contrasts)
  P <- stats::fitted(m)
  if (is.null(dim(P)) || ncol(P) == 1L) P <- cbind(1 - as.numeric(P), as.numeric(P))
  y <- stats::model.response(m$model)
  S <- do.call(cbind, lapply(seq.int(2, length(m$lev)), function(j)
    X * (as.numeric(y == m$lev[j]) - P[, j])))
  colnames(S) <- names(.sta_mnl_vector(m))
  S
}

.sta_iia_multinom <- function(model, methods, omit, seed, split, correction, control) {
  .sta_require("nnet")
  if (is.null(model$model)) stop("Fit multinom with model=TRUE.", call. = FALSE)
  if (!is.null(model$call$mask)) stop("IIA refitting does not support custom coefficient masks.", call. = FALSE)
  if (!is.null(model$decay) && model$decay != 0) stop("IIA refitting does not support penalized multinom.", call. = FALSE)
  weights <- stats::model.weights(model$model)
  if ((!is.null(weights) && any(weights != 1)) || !is.null(stats::model.offset(model$model)))
    stop("Multinomial IIA currently requires unweighted, offset-free models.", call. = FALSE)
  y <- stats::model.response(model$model)
  if (!is.factor(y) || is.ordered(y) || length(model$lev) < 3L)
    stop("IIA requires a nominal factor response with at least three categories.", call. = FALSE)
  if (!is.null(model$convergence) && model$convergence != 0) stop("The original multinom fit did not converge.", call. = FALSE)
  X <- stats::model.matrix(model$terms, model$model, contrasts.arg = model$contrasts)
  if (any(!is.finite(X)) || qr(X)$rank < ncol(X)) stop("IIA requires a finite, full-rank design matrix.", call. = FALSE)
  levels <- model$lev
  if (is.null(omit)) omit <- levels
  if (!is.character(omit) || !length(omit) || anyDuplicated(omit) || any(!omit %in% levels))
    stop("omit must contain distinct response category names.", call. = FALSE)
  defaults <- list(maxit = 1000L, abstol = 1e-12, reltol = 1e-12, MaxNWts = 10000L)
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control)) || any(!names(control) %in% names(defaults)))))
    stop("control may contain maxit, abstol, reltol and MaxNWts.", call. = FALSE)
  defaults[names(control)] <- control
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) previous <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed) assign(".Random.seed", previous, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  # Refit the original encoded design: factors, interactions and transformations
  # retain their exact fitted coding, without re-evaluating a raw-data formula.
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  fit <- function(rows, base) {
    yy <- droplevels(y[rows])
    expected <- unique(as.character(y[rows]))
    if (length(expected) < 2L || !base %in% expected) stop("A refit has fewer than two categories or lacks its base outcome.", call. = FALSE)
    yy <- stats::relevel(yy, base)
    d <- data.frame(.response = yy, X[rows, , drop = FALSE], check.names = FALSE)
    if (qr(X[rows, , drop = FALSE])$rank < ncol(X)) stop("A refit design is rank deficient; choose another partition or inspect rare levels.", call. = FALSE)
    m <- do.call(nnet::multinom, c(list(formula = .response ~ . - 1, data = d,
      Hess = TRUE, model = TRUE, trace = FALSE, decay = 0), defaults))
    if (m$convergence != 0L || any(!is.finite(stats::coef(m)))) stop("An IIA refit failed to converge.", call. = FALSE)
    m
  }
  n <- length(y); rows <- seq_len(n)
  if ("small_hsiao" %in% methods) {
    if (is.null(split)) {
      if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed))
        stop("seed must be a nonnegative finite integer.", call. = FALSE)
      set.seed(as.integer(seed))
      first <- sample(rows, floor(n/2))
      split <- rep(2L, n); split[first] <- 1L
    }
    if (length(split) != n || anyNA(split) || !setequal(unique(split), c(1,2)))
      stop("split must contain 1 or 2 for every fitted observation, with both halves nonempty.", call. = FALSE)
    first <- which(split == 1); second <- which(split == 2)
    for (r in list(first, second)) if (!setequal(as.character(y[r]), levels))
      stop("Both Small-Hsiao halves must contain all response categories; provide another split.", call. = FALSE)
  }
  output <- list()
  for (omitted in omit) {
    retained <- setdiff(levels, omitted); base <- retained[1]
    full <- fit(rows, base)
    keep <- as.character(y) != omitted
    reduced <- fit(rows[keep], base)
    nm <- names(.sta_mnl_vector(reduced))
    d <- .sta_mnl_vector(reduced) - .sta_mnl_vector(full)[nm]
    Vf <- .sta_mnl_vcov(full); Vr <- .sta_mnl_vcov(reduced)
    if ("hausman" %in% methods) output[[length(output)+1L]] <-
      .sta_iia_row("Hausman-McFadden", omitted, d, Vr - Vf[nm, nm, drop = FALSE])
    if ("suest" %in% methods) {
      If <- .sta_mnl_scores(full) %*% Vf
      Ir <- matrix(0, n, length(nm), dimnames = list(NULL, nm))
      Ir[keep, ] <- .sta_mnl_scores(reduced) %*% Vr
      D <- Ir - If[, nm, drop = FALSE]
      V <- crossprod(D) * if (correction) n/(n-1) else 1
      output[[length(output)+1L]] <- .sta_iia_row("Hausman (SUEST)", omitted, d, V)
    }
    if ("small_hsiao" %in% methods) {
      m1 <- fit(first, base); m2 <- fit(second, base)
      B <- stats::coef(m1)/sqrt(2) + (1 - 1/sqrt(2))*stats::coef(m2)
      restricted_rows <- second[as.character(y[second]) != omitted]
      rest <- fit(restricted_rows, base)
      B <- rbind(stats::setNames(rep(0, ncol(X)), colnames(X)), B)
      rownames(B) <- m2$lev
      eta <- X[restricted_rows, , drop = FALSE] %*% t(B[retained, , drop = FALSE])
      eta <- eta - apply(eta, 1, max)
      P <- exp(eta)/rowSums(exp(eta))
      ix <- cbind(seq_along(restricted_rows), match(as.character(y[restricted_rows]), retained))
      q <- 2 * (as.numeric(stats::logLik(rest)) - sum(log(P[ix])))
      k <- length(.sta_mnl_vector(rest))
      output[[length(output)+1L]] <- data.frame(Test = "Small-Hsiao", Omitted = omitted,
        statistic = q, df = k, p.value = if (q < 0) NA_real_ else stats::pchisq(q, k, lower.tail = FALSE),
        Status = if (q < 0) "Negative statistic: inspect fit" else "OK", stringsAsFactors = FALSE)
    }
  }
  list(table = do.call(rbind, output), split = if ("small_hsiao" %in% methods) split else NULL,
       notes = c("H0: independence of irrelevant alternatives (IIA). Failure to reject does not establish IIA.",
                 if ("small_hsiao" %in% methods) "Small-Hsiao depends on its partition; equal seed numbers in R and Stata do not imply equal partitions.",
                 if ("suest" %in% methods) paste0("SUEST uses the joint influence-function covariance; N/(N-1) correction ", if (correction) "enabled." else "disabled.")))
}

#' Tests of Independence of Irrelevant Alternatives
#' @param model A multinom, ordinary/nested mlogit or clogit fit.
#' @param method auto (default), hausman, suest, small_hsiao or lr. A vector
#' requests several tests. auto runs all three for multinom and Hausman for
#' ordinary mlogit/clogit. lr compares nested mlogit with its ordinary restriction.
#' @param omit Alternative/category names to remove in turn. Default: all for
#' multinom; all except the base alternative for ordinary mlogit.
#' @param reduced_model An explicitly fitted restricted model. Required for
#' clogit Hausman and nested-m logit LR. For Hausman, use the same specification,
#' with a restricted alternative set and its noninformative cases removed.
#' @param terms Exact coefficient names to compare in pairwise Hausman tests.
#' Default compares all common coefficients. Not used for multinom refits or LR.
#' @param order reduced-full (conventional Hausman covariance difference) or
#' full-reduced, to reproduce explicitly reversed examples. Applies to pairwise
#' clogit/mlogit Hausman, not multinom's conventional procedure.
#' @param seed Seed for Small-Hsiao. The caller's RNG state is restored.
#' @param split Optional vector of 1/2 labels aligned to the fitted observations.
#' Returned in result$split to reproduce the partition exactly.
#' @param correction Apply N/(N-1) to the joint SUEST covariance.
#' @param control Multinom refit controls: maxit, abstol, reltol and MaxNWts.
#' @param digits Display decimal places.
#' @return Invisibly, a statator_iia result with table, notes and the SH partition.
#' @details Multinom supports arbitrary full-rank encoded regressors, including
#' factors and interactions, but requires an unweighted, unpenalized, offset-free
#' nominal response with at least three categories and a retained model frame.
#' Refits with missing categories, singular designs or nonconvergence fail explicitly.
#' For ordinary mlogit, only Hausman is implemented; SUEST and Small-Hsiao are not
#' silently substituted by another test. Negative/indefinite Hausman covariance
#' differences have no interpretable p-value. LR restrictions and sample equality
#' are checked by ext_summary. This is a diagnostic, not a proof of IIA.
#' @examples
#' if (requireNamespace("nnet", quietly=TRUE)) {
#'   set.seed(10)
#'   d <- data.frame(x=rnorm(300), choice=factor(sample(letters[1:3],300,TRUE)))
#'   m <- nnet::multinom(choice~x, d, trace=FALSE, Hess=TRUE, model=TRUE)
#'   iia_test(m, method="hausman")
#' }
#' @export
iia_test <- function(model, method = "auto", omit = NULL, reduced_model = NULL,
                     terms = NULL, order = c("reduced-full", "full-reduced"),
                     seed = 12345, split = NULL, correction = TRUE,
                     control = list(), digits = 3) {
  order <- match.arg(order)
  .sta_flag(correction, "correction")
  if (inherits(model, "ext_summary")) model <- model$model
  if (inherits(reduced_model, "ext_summary")) reduced_model <- reduced_model$model
  if (inherits(model,"mlogit")) .sta_require("mlogit")
  if (inherits(model,"clogit")) .sta_require("survival")
  if (!inherits(model, c("multinom", "mlogit", "clogit"))) stop("IIA tests require a multinom, mlogit or clogit model.", call. = FALSE)
  if (!is.character(method) || !length(method) || anyNA(method) || anyDuplicated(method) ||
      any(!method %in% c("auto", "hausman", "suest", "small_hsiao", "lr"))) stop("Unknown or duplicate IIA method.", call. = FALSE)
  if ("auto" %in% method && length(method) > 1L) stop("auto cannot be combined with other methods.", call. = FALSE)
  if (identical(method, "auto")) method <- if (inherits(model, "multinom")) c("hausman", "suest", "small_hsiao") else "hausman"
  if (inherits(model,"mlogit") && !is.null(model$nests) && !identical(method,"lr"))
    stop("For nested mlogit use method='lr' with an ordinary reduced_model.", call. = FALSE)
  if (!"small_hsiao" %in% method && !is.null(split)) stop("split applies to Small-Hsiao only.", call. = FALSE)
  if (!inherits(model, "multinom") && length(control)) stop("control applies to multinom refits only.", call. = FALSE)
  if (inherits(model, "multinom")) {
    if ("lr" %in% method || !is.null(reduced_model) || !is.null(terms) || order != "reduced-full")
      stop("Multinom uses automatic conventional refits; LR/reduced_model/terms/reversed order are not supported.", call. = FALSE)
    result <- .sta_iia_multinom(model, method, omit, seed, split, correction, control)
  } else if (identical(method, "lr")) {
    if (!is.null(omit) || !is.null(terms) || order != "reduced-full") stop("omit/terms/reversed order do not apply to the nested LR test.", call. = FALSE)
    if (!inherits(model, "mlogit") || is.null(reduced_model)) stop("LR requires a nested mlogit and an ordinary reduced_model.", call. = FALSE)
    utils::capture.output(s <- ext_summary(model, iia_model = reduced_model))
    g <- s$iia_test
    result <- list(table = data.frame(Test = "LR: nested vs ordinary logit", Omitted = "None",
      statistic = g$statistic, df = g$df, p.value = g$p.value, Status = "OK"),
      notes = "H0: all nest dissimilarities equal one. Inspect convergence and admissible dissimilarity values separately.", split = NULL)
  } else {
    if (!identical(method, "hausman")) stop("This model supports Hausman only (or nested-vs-ordinary LR for mlogit).", call. = FALSE)
    if (!is.null(stats::model.weights(model$model)) || !is.null(stats::model.offset(model$model)))
      stop("Choice-model Hausman currently requires unweighted, offset-free fits.", call. = FALSE)
    if (!is.null(model$naive.var)) stop("Hausman requires model-based, not robust, fitted covariance.", call. = FALSE)
    compare <- function(full, restricted, omitted) {
      nm <- if (is.null(terms)) intersect(names(stats::coef(full)), names(stats::coef(restricted))) else terms
      if (!is.character(nm) || !length(nm) || anyDuplicated(nm) ||
          any(!nm %in% names(stats::coef(full))) || any(!nm %in% names(stats::coef(restricted))))
        stop("terms must identify common estimable coefficients.", call. = FALSE)
      d <- stats::coef(restricted)[nm] - stats::coef(full)[nm]
      V <- stats::vcov(restricted)[nm, nm, drop = FALSE] - stats::vcov(full)[nm, nm, drop = FALSE]
      if (order == "full-reduced") {d <- -d; V <- -V}
      .sta_iia_row("Hausman-McFadden", omitted, d, V)
    }
    if (!is.null(reduced_model)) {
      if (!inherits(reduced_model, class(model)[1])) stop("Full and restricted models must have the same fitted-model class.", call. = FALSE)
      if (!is.null(omit)) stop("Use omit for automatic mlogit refits OR reduced_model, not both.", call. = FALSE)
      if (!identical(attr(stats::terms(model), "term.labels"), attr(stats::terms(reduced_model), "term.labels")))
        stop("Hausman IIA requires the same regressors in full and restricted models.", call. = FALSE)
      result <- list(table = compare(model, reduced_model, "Supplied restriction"), split = NULL)
    } else {
      if (!inherits(model, "mlogit") || !is.null(model$nests)) stop("Supply reduced_model for clogit; automatic omission supports ordinary mlogit only.", call. = FALSE)
      info <- .sta_choice_data(model, NULL)
      model <- .sta_mlogit_call(model)
      reflevel <- .ext_mlogit(model)$base
      if (!length(reflevel)) reflevel <- info$alternatives[1]
      if (is.null(omit)) omit <- setdiff(info$alternatives, reflevel)
      if (!is.character(omit) || !length(omit) || anyDuplicated(omit) || any(!omit %in% info$alternatives)) stop("Unknown or duplicate omitted alternatives.", call. = FALSE)
      output <- lapply(omit, function(omitted) {
        retained <- setdiff(info$alternatives, omitted)
        if (length(retained) < 2L) stop("At least two alternatives must remain.", call. = FALSE)
        base <- if (reflevel %in% retained) reflevel else retained[1]
        full <- stats::update(model, data = model$model, reflevel = base)
        restricted <- stats::update(model, data = model$model, reflevel = base, alt.subset = retained)
        compare(full, restricted, omitted)
      })
      result <- list(table = do.call(rbind, output), split = NULL)
    }
    result$notes <- c("H0: independence of irrelevant alternatives. Non-rejection does not establish IIA.",
                       paste("Covariance difference order:", order),
                       "Supplied restrictions must remove alternatives, with the same utility specification; valid model pairing is the caller's responsibility.")
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0 || digits > 12 || digits != as.integer(digits)) stop("digits must be an integer from 0 to 12.", call. = FALSE)
  result$digits <- as.integer(digits)
  result$model_type <- class(model)[1]
  class(result) <- c("statator_iia", "statator_table")
  print(result)
  invisible(result)
}

#' @rdname iia_test
#' @param x An iia_test result.
#' @param ... Reserved for printing.
#' @export
print.statator_iia <- function(x, ...) {
  tab <- x$table
  display <- data.frame(Test = tab$Test, Omitted = tab$Omitted,
    Chi2 = ifelse(is.na(tab$statistic), ".", formatC(tab$statistic, format="f", digits=x$digits)),
    df = tab$df, `P>chi2` = ifelse(is.na(tab$p.value), ".", formatC(tab$p.value, format="f", digits=4)),
    Status = tab$Status, check.names=FALSE)
  cat("\nTests of independence of irrelevant alternatives (IIA)\n\n")
  cat("Model: ", x$model_type, "\n\n", sep="")
  display[] <- lapply(display, as.character)
  widths <- vapply(seq_along(display), function(j) max(nchar(c(names(display)[j], display[[j]]))), integer(1))
  line <- function(v) paste(vapply(seq_along(v), function(j) sprintf("%*s", widths[j], v[j]), character(1)), collapse="   ")
  rule <- strrep("-", sum(widths) + 3L*(length(widths)-1L))
  cat(rule, "\n", line(names(display)), "\n", rule, "\n", sep="")
  for (i in seq_len(nrow(display))) {
    if (i > 1L && display$Omitted[i] != display$Omitted[i-1L]) cat("\n")
    cat(line(unlist(display[i, ], use.names=FALSE)), "\n", sep="")
  }
  cat(rule, "\n", sep="")
  cat("\n")
  for (note in x$notes) cat("Note: ", note, "\n", sep="")
  cat("\n")
  invisible(x)
}
