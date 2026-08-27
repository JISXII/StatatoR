#' Box-Cox Regression
#'
#' Estimates Box-Cox regression models using maximum likelihood.
#'
#' @param formula Model formula.
#' @param data Data frame containing the variables.
#' @param model Box-Cox model specification.
#' @param notrans Variables that should not be transformed.
#' @param level Confidence level.
#' @param digits Number of decimal places.
#'
#' @return A Box-Cox regression result.
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
  # 0. VALIDACIONES
  # ============================================================

  model <- match.arg(model)

  if (length(level) != 1 || !is.numeric(level) ||
      level <= 0 || level >= 100) {
    stop("'level' debe estar entre 0 y 100.")
  }

  if (length(digits) != 1 || !is.numeric(digits) ||
      digits < 0 || digits != floor(digits)) {
    stop("'digits' debe ser un entero >= 0.")
  }

  # ============================================================
  # 1. MODEL FRAME
  # ============================================================

  mf <- model.frame(
    formula = formula,
    data = data,
    na.action = na.omit
  )

  y <- model.response(mf)

  if (!is.numeric(y)) {
    stop("La variable dependiente debe ser numérica.")
  }

  if (any(!is.finite(y))) {
    stop("La variable dependiente contiene valores no finitos.")
  }

  # ============================================================
  # 2. TERMINOS DEL MODELO
  # ============================================================

  trms <- terms(formula)

  rhs_labels <- attr(trms, "term.labels")

  X <- model.matrix(formula, mf)

  if (noconstant && "(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  assign_x <- attr(model.matrix(formula, mf), "assign")

  if (noconstant) {
    assign_x <- assign_x[
      colnames(X) != "(Intercept)"
    ]
  }

  # ============================================================
  # 3. IDENTIFICAR VARIABLES NO TRANSFORMADAS
  # ============================================================

  if (is.null(notrans)) {
    notrans <- character(0)
  }

  if (!is.character(notrans)) {
    stop("'notrans' debe ser un vector de nombres.")
  }

  all_vars <- all.vars(delete.response(trms))

  if (length(notrans) > 0 &&
      !all(notrans %in% all_vars)) {

    bad <- notrans[!notrans %in% all_vars]

    stop(
      paste(
        "Las siguientes variables de 'notrans' no están en el modelo:",
        paste(bad, collapse = ", ")
      )
    )
  }

  # ============================================================
  # 4. DETERMINAR COLUMNAS TRANSFORMABLES
  # ============================================================

  column_names <- colnames(X)

  transform_col <- rep(FALSE, ncol(X))

  if (ncol(X) > 0) {

    for (j in seq_along(column_names)) {

      nm <- column_names[j]

      if (nm == "(Intercept)") {
        transform_col[j] <- FALSE
        next
      }

      # detectar variable base
      base_name <- nm

      # intenta identificar nombres generados por factores
      matched <- all_vars[
        sapply(
          all_vars,
          function(v)
            nm == v ||
            startsWith(nm, v)
        )
      ]

      if (length(matched) == 1) {
        base_name <- matched
      }

      transform_col[j] <-
        !(base_name %in% notrans)
    }
  }

  # En lhsonly ninguna X se transforma
  if (model == "lhsonly") {
    transform_col[] <- FALSE
  }

  # ============================================================
  # 5. MATRICES
  # ============================================================

  Xtrans <- X[, transform_col, drop = FALSE]

  Xnotrans <- X[, !transform_col, drop = FALSE]

  # ============================================================
  # 6. POSITIVIDAD
  # ============================================================

  # Y se transforma en lhsonly, lambda y theta
  transform_y <- model %in% c(
    "lhsonly",
    "lambda",
    "theta"
  )

  if (transform_y && any(y <= 0)) {
    stop(
      "La variable dependiente debe ser estrictamente positiva ",
      "porque será transformada."
    )
  }

  # RHS se transforma en rhsonly, lambda y theta
  transform_x <- model %in% c(
    "rhsonly",
    "lambda",
    "theta"
  )

  if (transform_x && ncol(Xtrans) > 0) {

    if (any(Xtrans <= 0)) {

      bad_cols <- colnames(Xtrans)[
        apply(
          Xtrans,
          2,
          function(z) any(z <= 0)
        )
      ]

      stop(
        paste0(
          "Las siguientes variables contienen valores <= 0 ",
          "y serían transformadas: ",
          paste(bad_cols, collapse = ", "),
          ". Usa notrans = ... si corresponde."
        )
      )
    }
  }

  # ============================================================
  # 7. BOX-COX
  # ============================================================

  bc <- function(x, lambda) {

    lx <- log(x)

    if (abs(lambda) < 1e-8) {
      return(lx)
    }

    expm1(lambda * lx) / lambda
  }

  # ============================================================
  # 8. CONSTRUIR MODELO
  # ============================================================

  construct_model <- function(par) {

    if (model == "lhsonly") {

      theta <- par[1]

      lambda <- NA_real_

      yt <- bc(y, theta)

      W <- X

    } else if (model == "rhsonly") {

      lambda <- par[1]

      theta <- NA_real_

      yt <- y

      if (ncol(Xtrans) > 0) {

        Xt <- apply(
          Xtrans,
          2,
          function(z)
            bc(z, lambda)
        )

        Xt <- as.matrix(Xt)

        colnames(Xt) <- colnames(Xtrans)

      } else {
        Xt <- NULL
      }

      W <- cbind(
        if (!noconstant) rep(1, length(y)),
        Xt,
        Xnotrans
      )

      if (!noconstant) {
        colnames(W)[1] <- "(Intercept)"
      }

    } else if (model == "lambda") {

      lambda <- par[1]

      theta <- lambda

      yt <- bc(y, lambda)

      if (ncol(Xtrans) > 0) {

        Xt <- apply(
          Xtrans,
          2,
          function(z)
            bc(z, lambda)
        )

        Xt <- as.matrix(Xt)

        colnames(Xt) <- colnames(Xtrans)

      } else {
        Xt <- NULL
      }

      W <- cbind(
        if (!noconstant) rep(1, length(y)),
        Xt,
        Xnotrans
      )

      if (!noconstant) {
        colnames(W)[1] <- "(Intercept)"
      }

    } else {

      # theta model
      lambda <- par[1]
      theta <- par[2]

      yt <- bc(y, theta)

      if (ncol(Xtrans) > 0) {

        Xt <- apply(
          Xtrans,
          2,
          function(z)
            bc(z, lambda)
        )

        Xt <- as.matrix(Xt)

        colnames(Xt) <- colnames(Xtrans)

      } else {
        Xt <- NULL
      }

      W <- cbind(
        if (!noconstant) rep(1, length(y)),
        Xt,
        Xnotrans
      )

      if (!noconstant) {
        colnames(W)[1] <- "(Intercept)"
      }
    }

    list(
      y = yt,
      X = W,
      lambda = lambda,
      theta = theta
    )
  }

  # ============================================================
  # 9. LOG-LIKELIHOOD
  # ============================================================

  loglik <- function(par,
                     print_iter = FALSE,
                     prefix = "") {

    mod <- construct_model(par)

    fit <- tryCatch(
      lm.fit(mod$X, mod$y),
      error = function(e) NULL
    )

    if (is.null(fit))
      return(-Inf)

    res <- fit$residuals

    if (any(!is.finite(res)))
      return(-Inf)

    n <- length(y)

    k <- ncol(mod$X)

    if (n <= k)
      return(-Inf)

    SSR <- sum(res^2)

    if (!is.finite(SSR) || SSR <= 0)
      return(-Inf)

    sigma2 <- SSR / n

    ll <- -0.5 * n *
      (
        log(2 * pi) +
          1 +
          log(sigma2)
      )

    # ==========================================================
    # JACOBIAN
    # ==========================================================

    if (model == "lhsonly") {

      ll <- ll +
        (par[1] - 1) * sum(log(y))

    }

    if (model == "lambda") {

      ll <- ll +
        (par[1] - 1) * sum(log(y))

    }

    if (model == "theta") {

      ll <- ll +
        (par[2] - 1) * sum(log(y))

    }

    ll
  }

  # ============================================================
  # 10. COMPARISON MODEL
  # ============================================================

  # Restricted model: lambda/theta = 1
  #
  # Esto corresponde al modelo lineal original.
  # Se usa como comparación para LR chi2.

  comparison_par <-
    if (model == "theta") {
      c(1, 1)
    } else {
      1
    }

  ll_comparison <- loglik(
    comparison_par
  )

  # ============================================================
  # 11. STARTING VALUES
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
      c(1, 0.5),
      c(-0.5, 0),
      c(0, -0.5)
    )

  } else {

    starts <- matrix(
      c(
        -2,
        -1,
        -0.5,
        0,
        0.5,
        1,
        2
      ),
      ncol = 1
    )
  }

  # ============================================================
  # 12. OPTIMIZACIÓN CON TRACE
  # ============================================================

  optimize_model <- function(start) {

    iter <- 0

    objective <- function(p) {

      iter <<- iter + 1

      ll <- loglik(p)

      if (trace) {

        cat(
          "Iteration",
          iter - 1,
          ":  Log likelihood =",
          formatC(
            ll,
            digits = 8,
            format = "f"
          ),
          "\n"
        )
      }

      -ll
    }

    tryCatch(
      optim(
        par = start,
        fn = objective,
        method = "BFGS",
        hessian = TRUE,
        control = list(
          maxit = 5000,
          reltol = 1e-10
        )
      ),
      error = function(e) NULL
    )
  }

  # ============================================================
  # 13. ESTIMAR FULL MODEL
  # ============================================================

  if (trace) {
    cat("\n")
    cat("Fitting comparison model\n\n")
  }

  # mostrar comparison likelihood como punto inicial
  if (trace) {

    cat(
      "Iteration 0:  Log likelihood =",
      formatC(
        ll_comparison,
        digits = 8,
        format = "f"
      ),
      "\n"
    )
  }

  if (trace) {
    cat("\n")
    cat("Fitting full model\n\n")
  }

  fits <- lapply(
    seq_len(nrow(starts)),
    function(i) {

      # Para evitar repetir una salida enorme para cada starting value,
      # solo se usa trace completo en el primer intento.
      if (i > 1)
        trace_old <- trace

      optimize_model(starts[i, ])
    }
  )

  valid <- sapply(
    fits,
    function(z)
      !is.null(z) &&
      is.finite(z$value)
  )

  if (!any(valid)) {
    stop(
      "La optimización no pudo encontrar una solución válida."
    )
  }

  fits <- fits[valid]

  best <- fits[[
    which.min(
      sapply(
        fits,
        function(z) z$value
      )
    )
  ]]

  parhat <- best$par

  ll_full <- -best$value

  # ============================================================
  # 14. PARAMETROS DE TRANSFORMACION
  # ============================================================

  if (model == "lhsonly") {

    pnames <- "theta"

  } else if (model == "rhsonly") {

    pnames <- "lambda"

  } else if (model == "lambda") {

    pnames <- "lambda"

  } else {

    pnames <- c(
      "lambda",
      "theta"
    )
  }

  names(parhat) <- pnames

  # ============================================================
  # 15. VAR-COV DE PARAMETROS DE TRANSFORMACION
  # ============================================================

  H <- best$hessian

  V <- tryCatch(
    solve(H),
    error = function(e) NULL
  )

  if (is.null(V) ||
      any(!is.finite(V))) {

    V <- tryCatch(
      solve(
        optimHess(
          parhat,
          function(p)
            -loglik(p)
        )
      ),
      error = function(e) NULL
    )
  }

  if (!is.null(V)) {

    se_transform <- sqrt(
      pmax(
        diag(V),
        0
      )
    )

  } else {

    se_transform <-
      rep(
        NA_real_,
        length(parhat)
      )
  }

  # ============================================================
  # 16. Z, P Y IC
  # ============================================================

  z_transform <-
    parhat / se_transform

  p_transform <-
    2 * pnorm(
      abs(z_transform),
      lower.tail = FALSE
    )

  alpha <-
    1 - level / 100

  zcrit <-
    qnorm(
      1 - alpha / 2
    )

  lower <-
    parhat -
    zcrit * se_transform

  upper <-
    parhat +
    zcrit * se_transform

  transformation_table <-
    data.frame(
      Coefficient = parhat,
      `Std. err.` = se_transform,
      z = z_transform,
      `P>|z|` = p_transform,
      Lower = lower,
      Upper = upper,
      row.names = paste0(
        "/",
        names(parhat)
      ),
      check.names = FALSE
    )

  # ============================================================
  # 17. SCALE-VARIANT PARAMETERS
  # ============================================================

  mod_full <-
    construct_model(parhat)

  fit_full <-
    lm.fit(
      mod_full$X,
      mod_full$y
    )

  beta <- fit_full$coefficients

  residuals <- fit_full$residuals

  n <- length(y)

  sigma <-
    sqrt(
      sum(residuals^2) / n
    )

  beta_names <-
    colnames(mod_full$X)

  names(beta) <- beta_names

  # ============================================================
  # 18. LR TESTS
  # ============================================================

  restricted_values <-

    if (model == "theta") {

      list(
        "-1" = c(-1, -1),
        "0" = c(0, 0),
        "1" = c(1, 1)
      )

    } else {

      list(
        "-1" = -1,
        "0" = 0,
        "1" = 1
      )
    }

  lr_results <- list()

  for (nm in names(restricted_values)) {

    par0 <-
      restricted_values[[nm]]

    ll0 <-
      loglik(par0)

    df_lr <-
      length(parhat)

    LR <-
      2 * (
        ll_full -
          ll0
      )

    LR <-
      max(
        LR,
        0
      )

    p_lr <-
      pchisq(
        LR,
        df = df_lr,
        lower.tail = FALSE
      )

    lr_results[[nm]] <-
      c(
        Restricted =
          ll0,
        `LR statistic` =
          LR,
        df = df_lr,
        `Prob > chi2` =
          p_lr
      )
  }

  LR_table <-
    do.call(
      rbind,
      lr_results
    )

  rownames(LR_table) <-
    paste0(
      if (model == "theta")
        "theta = "
      else
        ifelse(
          model == "lhsonly",
          "theta = ",
          "lambda = "
        ),
      names(restricted_values)
    )

  # ============================================================
  # 19. LR DEL MODELO COMPLETO VS COMPARACIÓN
  # ============================================================

  LR_model <-
    2 * (
      ll_full -
        ll_comparison
    )

  df_model <-
    length(parhat)

  p_model <-
    pchisq(
      LR_model,
      df = df_model,
      lower.tail = FALSE
    )

  # ============================================================
  # 20. FORMATO
  # ============================================================

  fmt <- function(x) {

    ifelse(
      is.na(x),
      NA_character_,
      formatC(
        x,
        digits = digits,
        format = "f"
      )
    )
  }

  # ============================================================
  # 21. IMPRESION STATA-LIKE
  # ============================================================

  cat("\n")
  cat(
    "Number of obs   =",
    formatC(
      n,
      width = 10,
      format = "d"
    ),
    "\n"
  )

  cat(
    "LR chi2(",
    df_model,
    ")      =",
    formatC(
      LR_model,
      digits = 2,
      format = "f"
    ),
    "\n",
    sep = ""
  )

  cat(
    "Log likelihood =",
    fmt(ll_full),
    "\n"
  )

  cat(
    "Prob > chi2     =",
    ifelse(
      p_model < 0.001,
      "0.000",
      fmt(p_model)
    ),
    "\n\n"
  )

  # ============================================================
  # TABLA DE TRANSFORMACIÓN
  # ============================================================

  cat(
    "------------------------------------------------------------------------------\n"
  )

  dep_name <-
    as.character(
      formula[[2]]
    )

  cat(
    sprintf(
      "%12s | %12s %12s %8s %8s %12s %12s\n",
      dep_name,
      "Coefficient",
      "Std. err.",
      "z",
      "P>|z|",
      paste0(
        "[",
        level,
        "% conf. interval"
      ),
      ""
    )
  )

  cat(
    "-------------+----------------------------------------------------------------\n"
  )

  for (i in seq_len(nrow(transformation_table))) {

    cat(
      sprintf(
        "%12s | %12s %12s %8s %8s %12s %12s\n",
        rownames(
          transformation_table
        )[i],
        fmt(
          transformation_table$Coefficient[i]
        ),
        fmt(
          transformation_table$`Std. err.`[i]
        ),
        fmt(
          transformation_table$z[i]
        ),
        fmt(
          transformation_table$`P>|z|`[i]
        ),
        fmt(
          transformation_table$Lower[i]
        ),
        fmt(
          transformation_table$Upper[i]
        )
      )
    )
  }

  cat(
    "------------------------------------------------------------------------------\n\n"
  )

  # ============================================================
  # SCALE-VARIANT
  # ============================================================

  cat(
    "Estimates of scale-variant parameters\n"
  )

  cat(
    "----------------------------\n"
  )

  cat(
    "             | Coefficient\n"
  )

  cat(
    "-------------+--------------\n"
  )

  cat(
    "Notrans      |\n"
  )

  if (length(beta) > 0) {

    for (i in seq_along(beta)) {

      cat(
        sprintf(
          "%12s | %12s\n",
          names(beta)[i],
          fmt(beta[i])
        )
      )
    }
  }

  cat(
    "-------------+--------------\n"
  )

  cat(
    sprintf(
      "%12s | %12s\n",
      "/sigma",
      fmt(sigma)
    )
  )

  cat(
    "----------------------------\n\n"
  )

  # ============================================================
  # LR TESTS
  # ============================================================

  cat(
    "---------------------------------------------------------\n"
  )

  cat(
    "   Test         Restricted     LR statistic\n"
  )

  cat(
    "    H0:       log likelihood       chi2       Prob > chi2\n"
  )

  cat(
    "---------------------------------------------------------\n"
  )

  for (i in seq_len(nrow(LR_table))) {

    nm <-
      rownames(LR_table)[i]

    ll0 <-
      LR_table[i, "Restricted"]

    lr <-
      LR_table[i, "LR statistic"]

    pp <-
      LR_table[i, "Prob > chi2"]

    cat(
      sprintf(
        "%12s %16s %12s %15s\n",
        nm,
        fmt(ll0),
        fmt(lr),
        ifelse(
          pp < 0.001,
          "0.000",
          fmt(pp)
        )
      )
    )
  }

  cat(
    "---------------------------------------------------------\n"
  )

  # ============================================================
  # 22. OBJETO DE SALIDA
  # ============================================================

  result <- list(
    call = match.call(),
    model = model,
    n = n,
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
    iterations = best$counts["function"]
  )

  class(result) <-
    "boxcox_r"

  invisible(result)
}
