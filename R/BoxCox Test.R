#' Box-Cox Regression (replica robusta de `boxcox` de Stata)
#'
#' Estima modelos de regresion Box-Cox por maxima verosimilitud, replicando
#' exactamente las formulas de Stata (ver [R] boxcox, "Methods and formulas").
#'
#' @param formula Formula del modelo. Las variables del RHS que quieras
#'   transformar Y las que quieras dejar sin transformar (notrans) deben
#'   estar TODAS incluidas aca (a diferencia de Stata, donde notrans() agrega
#'   variables aunque no esten en indepvars).
#' @param data Data frame con las variables.
#' @param model Especificacion del modelo Box-Cox: "lhsonly", "rhsonly",
#'   "lambda" o "theta".
#' @param notrans Nombres de variables (tal como aparecen en `formula`, no en
#'   `colnames(model.matrix(...))`) que no deben transformarse.
#' @param level Nivel de confianza.
#' @param digits Numero de decimales para imprimir.
#' @param noconstant Si TRUE, omite el intercepto.
#' @param trace Si TRUE, imprime el log de iteraciones.
#'
#' @return Objeto de clase "boxcox_r" (invisible).
#' @export

boxcox <- function(formula,
                    data,
                    model = c("lhsonly", "rhsonly", "lambda", "theta"),
                    notrans = NULL,
                    level = 95,
                    digits = 4,
                    noconstant = FALSE,
                    trace = TRUE) {

  # ============================================================
  # 1. VALIDACIONES
  # ============================================================

  model <- match.arg(model)

  if (!is.numeric(level) || length(level) != 1 ||
      level <= 0 || level >= 100) {
    stop("'level' must be between 0 and 100.")
  }

  if (!is.numeric(digits) || length(digits) != 1 ||
      digits < 0 || digits != floor(digits)) {
    stop("'digits' must be a non-negative integer.")
  }

  # ============================================================
  # 2. MODEL FRAME
  # ============================================================

  mf <- model.frame(
    formula = formula,
    data = data,
    na.action = na.omit
  )

  y <- model.response(mf)

  if (!is.numeric(y)) {
    stop("The dependent variable must be numeric.")
  }

  if (model %in% c("lhsonly", "lambda", "theta")) {

    if (any(!is.finite(y))) {
      stop("The dependent variable contains non-finite values.")
    }

    if (any(y <= 0)) {
      stop(
        "The dependent variable must be strictly positive ",
        "for this Box-Cox model."
      )
    }
  }

  # ============================================================
  # 3. DESIGN MATRIX
  # ============================================================

  X <- model.matrix(formula, mf)

  # Guardamos el mapeo columna -> termino ANTES de tocar nada.
  # attr(X, "assign")[j] == 0  -> columna j es el intercepto
  # attr(X, "assign")[j] == k  -> columna j viene de term.labels(formula)[k]
  term_labels  <- attr(terms(formula), "term.labels")
  assign_vec   <- attr(X, "assign")

  if (noconstant && "(Intercept)" %in% colnames(X)) {
    keep       <- colnames(X) != "(Intercept)"
    X          <- X[, keep, drop = FALSE]
    assign_vec <- assign_vec[keep]
  }

  # ============================================================
  # 4. NOTRANS
  # ============================================================

  if (is.null(notrans)) {
    notrans <- character(0)
  }

  if (!is.character(notrans)) {
    stop("'notrans' must contain variable names.")
  }

  if (length(notrans) > 0 && !all(notrans %in% term_labels)) {
    faltantes <- notrans[!notrans %in% term_labels]
    stop(
      "'notrans' contains variables not present as terms in 'formula': ",
      paste(faltantes, collapse = ", "),
      ". Include them in 'formula' first (see documentation)."
    )
  }

  # ============================================================
  # 5. IDENTIFY TRANSFORMED RHS VARIABLES
  #    (mapeo EXACTO via attr(X,"assign"), no coincidencia de texto)
  # ============================================================

  transform_x <- model %in% c("rhsonly", "lambda", "theta")

  if (transform_x) {

    transform_cols <- vapply(
      assign_vec,
      function(a) {
        if (a == 0) return(FALSE)                 # intercepto: nunca
        !(term_labels[a] %in% notrans)             # transformar si NO esta en notrans
      },
      logical(1)
    )

  } else {

    transform_cols <- rep(FALSE, ncol(X))
  }

  # ============================================================
  # 6. CHECK POSITIVITY OF RHS
  # ============================================================

  if (transform_x && any(transform_cols)) {

    rhs_vars <- all.vars(delete.response(terms(formula)))

    for (v in rhs_vars) {

      if (v %in% notrans) next
      if (!v %in% names(mf)) next

      xv <- mf[[v]]

      if (is.numeric(xv) && any(xv <= 0)) {
        stop(
          paste0(
            "Variable '", v, "' contains values <= 0. ",
            "It must be strictly positive because it is transformed."
          )
        )
      }
    }
  }

  # ============================================================
  # 7. BOX-COX TRANSFORMATION (mismo umbral que Stata: 1e-10)
  # ============================================================

  BC_EPS <- 1e-10

  bc <- function(x, lambda) {
    lx <- log(x)
    if (abs(lambda) < BC_EPS) {
      return(lx)
    } else {
      return(expm1(lambda * lx) / lambda)
    }
  }

  # ============================================================
  # 8. CONSTRUCT MODEL
  # ============================================================

  construct_model <- function(par) {

    if (model == "lhsonly") {

      theta <- par[1]
      yt <- bc(y, theta)
      W <- X

      return(list(y = yt, X = W, lambda = NA_real_, theta = theta))
    }

    if (model == "rhsonly") {

      lambda <- par[1]
      yt <- y
      W <- X

      for (j in seq_len(ncol(W))) {
        if (transform_cols[j]) W[, j] <- bc(W[, j], lambda)
      }

      return(list(y = yt, X = W, lambda = lambda, theta = NA_real_))
    }

    if (model == "lambda") {

      lambda <- par[1]
      yt <- bc(y, lambda)
      W <- X

      for (j in seq_len(ncol(W))) {
        if (transform_cols[j]) W[, j] <- bc(W[, j], lambda)
      }

      return(list(y = yt, X = W, lambda = lambda, theta = lambda))
    }

    if (model == "theta") {

      lambda <- par[1]
      theta  <- par[2]
      yt <- bc(y, theta)
      W <- X

      for (j in seq_len(ncol(W))) {
        if (transform_cols[j]) W[, j] <- bc(W[, j], lambda)
      }

      return(list(y = yt, X = W, lambda = lambda, theta = theta))
    }
  }

  # ============================================================
  # 9. CONCENTRATED LOG-LIKELIHOOD
  #    (identica a la de Stata: [R] boxcox, Methods and formulas)
  # ============================================================

  n <- length(y)
  sum_log_y <- sum(log(y))   # se recalcula una sola vez

  loglik <- function(par) {

    mod <- construct_model(par)

    fit <- tryCatch(
      lm.fit(mod$X, mod$y),
      error = function(e) NULL
    )

    if (is.null(fit)) return(-Inf)

    # Salvaguarda: si el diseño quedo colineal/deficiente en rango para
    # este valor de lambda/theta, la verosimilitud NO es comparable con
    # la de un ajuste de rango completo -> descartamos ese punto.
    if (fit$qr$rank < ncol(mod$X)) return(-Inf)

    residuals <- fit$residuals
    SSR <- sum(residuals^2)

    if (!is.finite(SSR) || SSR <= 0) return(-Inf)

    sigma2 <- SSR / n

    ll <- -0.5 * n * (log(2 * pi) + 1 + log(sigma2))

    if (model == "lhsonly") ll <- ll + (par[1] - 1) * sum_log_y
    if (model == "lambda")  ll <- ll + (par[1] - 1) * sum_log_y
    if (model == "theta")   ll <- ll + (par[2] - 1) * sum_log_y
    # rhsonly: sin termino Jacobiano (y no se transforma)

    ll
  }

  # ============================================================
  # 10. COMPARISON MODEL (lambda = theta = 1)
  # ============================================================

  comparison_par <- if (model == "theta") c(1, 1) else 1
  ll_comparison  <- loglik(comparison_par)

  # ============================================================
  # 11. STARTING VALUES
  #     Igual que Stata: arranca en 1. Ademas, para mayor robustez
  #     (superficies planas/tipo cresta, sobre todo en "theta"),
  #     probamos una grilla adicional y nos quedamos con el mejor.
  # ============================================================

  if (model == "theta") {
    starts <- rbind(
      c(1, 1),      # arranque de Stata
      c(0, 0),
      c(0.5, 0.5),
      c(-1, -1),
      c(1, 0),
      c(0, 1),
      c(0.5, 1),
      c(1, 0.5)
    )
  } else {
    starts <- matrix(
      c(1, -2, -1, -0.5, 0, 0.5, 2),   # 1 primero: arranque de Stata
      ncol = 1
    )
  }

  # ============================================================
  # 12. OPTIMIZACION
  #     BFGS (grilla de arranques) + refinamiento final con Nelder-Mead
  #     y tolerancias estrictas, para converger al mismo optimo que
  #     el Newton-Raphson de Stata.
  # ============================================================

  if (trace) {
    cat("\n")
    cat("Fitting comparison model\n\n")
    cat(
      "Iteration 0:  Log likelihood = ",
      formatC(ll_comparison, digits = 8, format = "f"),
      "\n", sep = ""
    )
    cat("\nFitting full model\n\n")
  }

  best <- NULL
  best_value <- Inf

  fit_one_start <- function(par0, verbose) {

    current_iter <- 0

    objective <- function(par) {
      value <- -loglik(par)
      current_iter <<- current_iter + 1
      if (verbose) {
        cat(
          "Iteration ", current_iter - 1, ":  Log likelihood = ",
          formatC(-value, digits = 8, format = "f"), "\n", sep = ""
        )
      }
      value
    }

    fit1 <- tryCatch(
      optim(
        par = par0, fn = objective, method = "BFGS",
        control = list(maxit = 5000, reltol = 1e-12)
      ),
      error = function(e) NULL
    )

    if (is.null(fit1) || !is.finite(fit1$value)) return(NULL)

    # Refinamiento: Nelder-Mead no usa gradiente numerico y suele
    # limar los ultimos decimales que BFGS deja sueltos en superficies
    # chatas (muy comun en el modelo "theta").
    fit2 <- tryCatch(
      optim(
        par = fit1$par, fn = objective, method = "Nelder-Mead",
        control = list(maxit = 5000, reltol = 1e-14)
      ),
      error = function(e) NULL
    )

    if (!is.null(fit2) && is.finite(fit2$value) && fit2$value <= fit1$value) {
      fit1 <- fit2
    }

    # Hessiano final siempre por diferencias finitas de alta precision,
    # independiente del metodo que gano (mas estable que hessian=TRUE
    # de optim con BFGS cuando el refinamiento fue Nelder-Mead).
    fit1
  }

  for (s in seq_len(nrow(starts))) {

    fit <- fit_one_start(starts[s, ], verbose = (s == 1 && trace))

    if (!is.null(fit) && fit$value < best_value) {
      best <- fit
      best_value <- fit$value
    }

    if (s == 1 && trace) cat("\n")
  }

  if (is.null(best)) {
    stop("Maximum likelihood optimization failed.")
  }

  parhat <- best$par
  ll_full <- -best$value

  # ============================================================
  # 13. PARAMETER NAMES
  # ============================================================

  if (model == "lhsonly") {
    names(parhat) <- "theta"
  } else if (model == "rhsonly") {
    names(parhat) <- "lambda"
  } else if (model == "lambda") {
    names(parhat) <- "lambda"
  } else {
    names(parhat) <- c("lambda", "theta")
  }

  # ============================================================
  # 14. VARIANCE-COVARIANCE MATRIX
  #     Hessiano numerico de alta precision (optimHess), consistente
  #     sin importar que optimizador gano el refinamiento.
  # ============================================================

  H <- tryCatch(
    optimHess(parhat, function(p) -loglik(p),
              control = list(ndeps = rep(1e-6, length(parhat)))),
    error = function(e) NULL
  )

  V <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL

  if (!is.null(V)) {
    se_transform <- sqrt(pmax(diag(V), 0))
  } else {
    se_transform <- rep(NA_real_, length(parhat))
  }

  # ============================================================
  # 15. CONFIDENCE INTERVAL
  # ============================================================

  alpha <- 1 - (level / 100)
  zcrit <- qnorm(1 - alpha / 2)

  z_transform <- parhat / se_transform
  p_transform <- 2 * pnorm(abs(z_transform), lower.tail = FALSE)

  lower <- parhat - zcrit * se_transform
  upper <- parhat + zcrit * se_transform

  transformation_table <- data.frame(
    Coefficient = parhat,
    `Std. err.` = se_transform,
    z = z_transform,
    `P>|z|` = p_transform,
    Lower = lower,
    Upper = upper,
    row.names = paste0("/", names(parhat)),
    check.names = FALSE
  )

  # ============================================================
  # 16. SCALE-VARIANT PARAMETERS
  # ============================================================

  mod_full <- construct_model(parhat)
  fit_full <- lm.fit(mod_full$X, mod_full$y)

  if (fit_full$qr$rank < ncol(mod_full$X)) {
    warning(
      "The design matrix is rank-deficient at the optimum ",
      "(collinearity after Box-Cox transformation). Coefficients ",
      "for the aliased columns are NA."
    )
  }

  beta  <- fit_full$coefficients
  sigma <- sqrt(sum(fit_full$residuals^2) / length(y))

  names(beta) <- colnames(mod_full$X)

  # ============================================================
  # 17. LR TESTS
  # ============================================================

  restricted_values <- if (model == "theta") {
    list("-1" = c(-1, -1), "0" = c(0, 0), "1" = c(1, 1))
  } else {
    list("-1" = -1, "0" = 0, "1" = 1)
  }

  LR_tests <- list()

  for (nm in names(restricted_values)) {

    ll0 <- loglik(restricted_values[[nm]])
    LR  <- max(0, 2 * (ll_full - ll0))
    df_lr <- length(parhat)
    p_lr  <- pchisq(LR, df = df_lr, lower.tail = FALSE)

    LR_tests[[nm]] <- c(
      Restricted = ll0,
      `LR statistic` = LR,
      df = df_lr,
      `Prob > chi2` = p_lr
    )
  }

  LR_table <- do.call(rbind, LR_tests)

  test_parameter <- if (model == "lhsonly" || model == "theta") "theta" else "lambda"

  rownames(LR_table) <- paste0(test_parameter, " = ", names(restricted_values))

  # ============================================================
  # 18. MODEL LR TEST
  # ============================================================

  LR_model <- max(0, 2 * (ll_full - ll_comparison))
  df_model <- length(parhat)
  p_model  <- pchisq(LR_model, df = df_model, lower.tail = FALSE)

  # ============================================================
  # 19. FORMATTING
  # ============================================================

  fmt <- function(x) {
    ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
  }

  fmt_p <- function(x) {
    ifelse(x < 0.001, "0.000", fmt(x))
  }

  # ============================================================
  # 20. HEADER
  # ============================================================

  cat("Number of obs   = ", nrow(mf), "\n", sep = "")
  cat("Log likelihood = ", fmt(ll_full), "\n", sep = "")

  # ============================================================
  # 21. TRANSFORMATION TABLE
  # ============================================================

  cat("------------------------------------------------------------------------------\n")

  dep_name <- all.vars(formula)[1]
  if (length(dep_name) != 1 || is.na(dep_name)) dep_name <- "y"

  cat(
    sprintf(
      "%12s | %12s %12s %8s %8s %12s %12s\n",
      dep_name, "Coefficient", "Std. err.", "z", "P>|z|",
      paste0("[", level, "% conf. interval"), ""
    )
  )

  cat("-------------+----------------------------------------------------------------\n")

  for (i in seq_len(nrow(transformation_table))) {
    cat(
      sprintf(
        "%12s | %12s %12s %8s %8s %12s %12s\n",
        rownames(transformation_table)[i],
        fmt(transformation_table$Coefficient[i]),
        fmt(transformation_table$`Std. err.`[i]),
        fmt(transformation_table$z[i]),
        fmt_p(transformation_table$`P>|z|`[i]),
        fmt(transformation_table$Lower[i]),
        fmt(transformation_table$Upper[i])
      )
    )
  }

  cat("------------------------------------------------------------------------------\n\n")

  # ============================================================
  # 22. SCALE-VARIANT PARAMETERS
  # ============================================================

  cat("Estimates of scale-variant parameters\n")
  cat("----------------------------\n")
  cat("             | Coefficient\n")
  cat("-------------+--------------\n")
  cat("Notrans      |\n")

  if (length(beta) > 0) {
    for (i in seq_along(beta)) {
      cat(sprintf("%12s | %12s\n", names(beta)[i], fmt(beta[i])))
    }
  }

  cat("-------------+--------------\n")
  cat(sprintf("%12s | %12s\n", "/sigma", fmt(sigma)))
  cat("----------------------------\n\n")

  # ============================================================
  # 23. LR TESTS
  # ============================================================

  cat("---------------------------------------------------------\n")
  cat("   Test         Restricted     LR statistic\n")
  cat("    H0:       log likelihood       chi2       Prob > chi2\n")
  cat("---------------------------------------------------------\n")

  for (i in seq_len(nrow(LR_table))) {
    cat(
      sprintf(
        "%12s %16s %12s %15s\n",
        rownames(LR_table)[i],
        fmt(LR_table[i, "Restricted"]),
        fmt(LR_table[i, "LR statistic"]),
        fmt_p(LR_table[i, "Prob > chi2"])
      )
    )
  }

  cat("---------------------------------------------------------\n")

  # ============================================================
  # 24. RETURN OBJECT
  # ============================================================

  result <- list(
    call = match.call(),
    model = model,
    n = nrow(mf),
    logLik = ll_full,
    comparison_logLik = ll_comparison,
    LR = LR_model,
    df = df_model,
    p_value = p_model,
    transformation = transformation_table,
    coefficients = beta,
    sigma = sigma,
    LR_tests = LR_table,
    vcov = V,
    hessian = H,
    convergence = best$convergence,
    rank_full = fit_full$qr$rank,
    ncol_full = ncol(mod_full$X)
  )

  class(result) <- "boxcox_r"

  invisible(result)
}
