#' Box-Cox Regression (replica robusta de `boxcox` de Stata)
#'
#' Estima modelos de regresion Box-Cox por maxima verosimilitud, replicando
#' las formulas y el comportamiento de Stata (ver [R] boxcox, "Methods and
#' formulas"), incluyendo el manejo de colinealidad estatica que Stata
#' reporta como "note: <var> omitted because of collinearity."
#'
#' @param formula Formula del modelo. Las variables que quieras transformar
#'   Y las que quieras dejar en notrans deben estar TODAS incluidas aca
#'   (a diferencia de Stata, donde notrans() puede agregar variables aunque
#'   no esten en indepvars).
#' @param data Data frame con las variables.
#' @param model "lhsonly", "rhsonly", "lambda" o "theta".
#' @param notrans Nombres de variables (tal como aparecen en `formula`) que
#'   no deben transformarse. Si alguna resulta colineal con el resto del
#'   bloque no-transformado (intercepto + otras notrans), se reclasifica
#'   automaticamente como transformada -igual que hace Stata- y se avisa
#'   por consola. Solo aplica a variables numericas continuas; si la
#'   colineal es un factor o el intercepto, se detiene con un error.
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

  # Mapeo exacto columna -> termino del formula (0 = intercepto).
  term_labels <- attr(terms(formula), "term.labels")
  assign_vec  <- attr(X, "assign")

  if (noconstant && "(Intercept)" %in% colnames(X)) {
    keep       <- colnames(X) != "(Intercept)"
    X          <- X[, keep, drop = FALSE]
    assign_vec <- assign_vec[keep]
  }

  # ============================================================
  # 4. NOTRANS: validacion de nombres
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

  transform_x <- model %in% c("rhsonly", "lambda", "theta")

  # ============================================================
  # 5. CHEQUEO ESTATICO DE COLINEALIDAD DEL BLOQUE NO TRANSFORMADO
  #    (equivalente al "note: x1 omitted because of collinearity."
  #    de Stata)
  #
  #    El bloque [intercepto + variables en notrans] nunca cambia
  #    durante la optimizacion (a diferencia de las columnas que se
  #    transforman con lambda/theta). Si ese bloque es de rango
  #    deficiente, el modelo no es identificable pase lo que pase con
  #    lambda/theta. Stata resuelve esto sacando la variable ofensora
  #    de notrans y dejando que SI se transforme - eso es justamente
  #    lo que reclasificamos aca.
  # ============================================================

  notrans_effective <- notrans

  if (transform_x && length(notrans_effective) > 0) {

    repeat {

      static_cols <- vapply(
        seq_len(ncol(X)),
        function(j) {
          a <- assign_vec[j]
          if (a == 0) return(TRUE)
          term_labels[a] %in% notrans_effective
        },
        logical(1)
      )

      static_idx <- which(static_cols)
      Z <- X[, static_idx, drop = FALSE]

      qrZ   <- qr(Z)
      rankZ <- qrZ$rank

      if (rankZ >= ncol(Z)) break   # sin colinealidad, listo

      redundant_local <- qrZ$pivot[(rankZ + 1):ncol(Z)]
      redundant_cols  <- static_idx[redundant_local]

      moved_any <- FALSE

      for (col in redundant_cols) {

        if (assign_vec[col] == 0) {
          stop(
            "The intercept is collinear with the untransformed block. ",
            "This model is not identifiable as specified; check your data."
          )
        }

        term_nm <- term_labels[assign_vec[col]]

        if (!(term_nm %in% notrans_effective)) next  # ya fue reclasificado antes

        # Solo se puede "salvar" reclasificando si el termino es numerico
        # continuo (no factor) y estrictamente positivo, porque va a
        # requerir Box-Cox.
        raw_var <- if (term_nm %in% names(mf)) mf[[term_nm]] else NULL

        if (is.null(raw_var) || !is.numeric(raw_var) || is.factor(raw_var)) {
          stop(
            "Variable '", term_nm, "' is collinear with the rest of the ",
            "untransformed block, and Stata would normally drop it from ",
            "notrans() and transform it instead. That fallback only works ",
            "for continuous numeric variables here (this one is categorical ",
            "or non-numeric), so please resolve the collinearity manually ",
            "(e.g. remove the variable or the redundant one)."
          )
        }

        if (any(raw_var <= 0)) {
          stop(
            "Variable '", term_nm, "' is collinear with the rest of the ",
            "untransformed block. Stata would drop it from notrans() and ",
            "transform it instead, but it contains values <= 0 so it ",
            "cannot receive a Box-Cox transform. Resolve the collinearity ",
            "manually."
          )
        }

        cat(
          "note: ", term_nm,
          " is collinear with the untransformed block; ",
          "reclassifying it from notrans to transformed (as Stata does).\n",
          sep = ""
        )

        notrans_effective <- setdiff(notrans_effective, term_nm)
        moved_any <- TRUE
      }

      if (!moved_any) {
        stop(
          "The untransformed block (intercept + notrans variables) is ",
          "rank-deficient and the collinearity could not be resolved ",
          "automatically. Check your data for exact linear dependencies."
        )
      }
    }
  }

  # ============================================================
  # 6. IDENTIFY TRANSFORMED RHS VARIABLES (definitivo, tras el
  #    chequeo de colinealidad)
  # ============================================================

  if (transform_x) {

    transform_cols <- vapply(
      seq_len(ncol(X)),
      function(j) {
        a <- assign_vec[j]
        if (a == 0) return(FALSE)
        !(term_labels[a] %in% notrans_effective)
      },
      logical(1)
    )

  } else {

    transform_cols <- rep(FALSE, ncol(X))
  }

  # ============================================================
  # 7. CHECK POSITIVITY OF TRANSFORMED RHS VARIABLES
  # ============================================================

  if (transform_x && any(transform_cols)) {

    rhs_vars <- all.vars(delete.response(terms(formula)))

    for (v in rhs_vars) {

      if (v %in% notrans_effective) next
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
  # 8. BOX-COX TRANSFORMATION (mismo umbral que Stata: 1e-10)
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
  # 9. CONSTRUCT MODEL
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
  # 10. CONCENTRATED LOG-LIKELIHOOD
  #     (identica a la de Stata: [R] boxcox, Methods and formulas)
  #
  #     IMPORTANTE: NO se penaliza con -Inf la deficiencia de rango.
  #     lm.fit() maneja columnas colineales via QR con pivoteo: los
  #     valores ajustados (y por lo tanto el SSR y la verosimilitud)
  #     siguen siendo correctos aunque coeficientes individuales no
  #     esten identificados. Forzar -Inf ahi (como se hacia antes)
  #     introduce un acantilado discontinuo en la superficie que
  #     rompe el gradiente numerico de BFGS - eso era exactamente lo
  #     que generaba las oscilaciones y el "-Inf" espurio del log.
  # ============================================================

  n <- length(y)
  sum_log_y <- sum(log(y))

  loglik <- function(par) {

    mod <- construct_model(par)

    fit <- tryCatch(
      lm.fit(mod$X, mod$y),
      error = function(e) NULL
    )

    if (is.null(fit)) return(-Inf)

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
  # 11. COMPARISON MODEL (lambda = theta = 1)
  # ============================================================

  comparison_par <- if (model == "theta") c(1, 1) else 1
  ll_comparison  <- loglik(comparison_par)

  # ============================================================
  # 12. STARTING VALUES (arranque de Stata primero: 1)
  # ============================================================

  if (model == "theta") {
    starts <- rbind(
      c(1, 1),
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
      c(1, -2, -1, -0.5, 0, 0.5, 2),
      ncol = 1
    )
  }

  # ============================================================
  # 13. OPTIMIZACION: BFGS (grilla de arranques) + refinamiento
  #     final con Nelder-Mead y tolerancias estrictas.
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
  # 14. PARAMETER NAMES
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
  # 15. VARIANCE-COVARIANCE MATRIX
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
  # 16. CONFIDENCE INTERVAL
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
  # 17. SCALE-VARIANT PARAMETERS
  # ============================================================

  mod_full <- construct_model(parhat)
  fit_full <- lm.fit(mod_full$X, mod_full$y)

  if (fit_full$qr$rank < ncol(mod_full$X)) {
    warning(
      "The design matrix is rank-deficient at the optimum. Coefficients ",
      "for the aliased columns are NA. This can happen at pathological ",
      "lambda/theta values even after the static collinearity fix; the ",
      "fitted values and log-likelihood remain valid."
    )
  }

  beta  <- fit_full$coefficients
  sigma <- sqrt(sum(fit_full$residuals^2) / length(y))

  names(beta) <- colnames(mod_full$X)

  # ============================================================
  # 18. LR TESTS
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
  # 19. MODEL LR TEST
  # ============================================================

  LR_model <- max(0, 2 * (ll_full - ll_comparison))
  df_model <- length(parhat)
  p_model  <- pchisq(LR_model, df = df_model, lower.tail = FALSE)

  # ============================================================
  # 20. FORMATTING
  # ============================================================

  fmt <- function(x) {
    ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
  }

  fmt_p <- function(x) {
    ifelse(x < 0.001, "0.000", fmt(x))
  }

  # ============================================================
  # 21. HEADER
  # ============================================================

  cat("Number of obs   = ", nrow(mf), "\n", sep = "")
  cat("Log likelihood = ", fmt(ll_full), "\n", sep = "")

  # ============================================================
  # 22. TRANSFORMATION TABLE
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
  # 23. SCALE-VARIANT PARAMETERS
  #     (separado en Notrans / Trans como hace Stata de verdad; la
  #     version anterior mostraba todo bajo "Notrans" siempre, lo
  #     cual deja de tener sentido apenas hay columnas Trans)
  # ============================================================

  cat("Estimates of scale-variant parameters\n")
  cat("----------------------------\n")
  cat("             | Coefficient\n")
  cat("-------------+--------------\n")

  is_notrans_col <- !transform_cols   # incluye intercepto
  is_trans_col   <- transform_cols

  if (any(is_notrans_col)) {
    cat("Notrans      |\n")
    for (i in which(is_notrans_col)) {
      cat(sprintf("%12s | %12s\n", names(beta)[i], fmt(beta[i])))
    }
  }

  if (any(is_trans_col)) {
    if (any(is_notrans_col)) cat("-------------+--------------\n")
    cat("Trans        |\n")
    for (i in which(is_trans_col)) {
      cat(sprintf("%12s | %12s\n", names(beta)[i], fmt(beta[i])))
    }
  }

  cat("-------------+--------------\n")
  cat(sprintf("%12s | %12s\n", "/sigma", fmt(sigma)))
  cat("----------------------------\n\n")

  # ============================================================
  # 24. LR TESTS
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
  # 25. RETURN OBJECT
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
    notrans_effective = notrans_effective,
    transform_cols = transform_cols,
    rank_full = fit_full$qr$rank,
    ncol_full = ncol(mod_full$X)
  )

  class(result) <- "boxcox_r"

  invisible(result)
}
