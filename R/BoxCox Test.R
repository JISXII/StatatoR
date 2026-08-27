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

  if (noconstant && "(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
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

  # ============================================================
  # 5. IDENTIFY TRANSFORMED RHS VARIABLES
  # ============================================================

  transform_x <- model %in% c(
    "rhsonly",
    "lambda",
    "theta"
  )

  if (transform_x && length(notrans) == 0) {

    transform_cols <- colnames(X) != "(Intercept)"

  } else if (transform_x) {

    transform_cols <- sapply(
      colnames(X),
      function(nm) {

        if (nm == "(Intercept)")
          return(FALSE)

        base_name <- nm

        matches <- notrans[
          sapply(
            notrans,
            function(v) {
              nm == v ||
                startsWith(nm, paste0(v))
            }
          )
        ]

        length(matches) == 0
      }
    )

  } else {

    transform_cols <- rep(
      FALSE,
      ncol(X)
    )
  }

  # ============================================================
  # 6. CHECK POSITIVITY OF RHS
  # ============================================================

  if (transform_x &&
      any(transform_cols)) {

    # Identify original variables in formula
    rhs_vars <- all.vars(
      delete.response(
        terms(formula)
      )
    )

    for (v in rhs_vars) {

      if (v %in% notrans)
        next

      if (!v %in% names(mf))
        next

      xv <- mf[[v]]

      if (is.numeric(xv) &&
          any(xv <= 0)) {

        stop(
          paste0(
            "Variable '",
            v,
            "' contains values <= 0. ",
            "It must be strictly positive because it is transformed."
          )
        )
      }
    }
  }

  # ============================================================
  # 7. BOX-COX TRANSFORMATION
  # ============================================================

  bc <- function(x, lambda) {

    lx <- log(x)

    if (abs(lambda) < 1e-8) {

      return(lx)

    } else {

      return(
        expm1(lambda * lx) / lambda
      )
    }
  }

  # ============================================================
  # 8. CONSTRUCT MODEL
  # ============================================================

  construct_model <- function(par) {

    # ----------------------------------------------------------
    # LHS ONLY
    # ----------------------------------------------------------

    if (model == "lhsonly") {

      theta <- par[1]

      yt <- bc(
        y,
        theta
      )

      W <- X

      return(
        list(
          y = yt,
          X = W,
          lambda = NA_real_,
          theta = theta
        )
      )
    }

    # ----------------------------------------------------------
    # RHS ONLY
    # ----------------------------------------------------------

    if (model == "rhsonly") {

      lambda <- par[1]

      yt <- y

      W <- X

      for (j in seq_len(ncol(W))) {

        if (transform_cols[j]) {

          W[, j] <- bc(
            W[, j],
            lambda
          )
        }
      }

      return(
        list(
          y = yt,
          X = W,
          lambda = lambda,
          theta = NA_real_
        )
      )
    }

    # ----------------------------------------------------------
    # SAME PARAMETER BOTH SIDES
    # ----------------------------------------------------------

    if (model == "lambda") {

      lambda <- par[1]

      yt <- bc(
        y,
        lambda
      )

      W <- X

      for (j in seq_len(ncol(W))) {

        if (transform_cols[j]) {

          W[, j] <- bc(
            W[, j],
            lambda
          )
        }
      }

      return(
        list(
          y = yt,
          X = W,
          lambda = lambda,
          theta = lambda
        )
      )
    }

    # ----------------------------------------------------------
    # DIFFERENT PARAMETERS
    # ----------------------------------------------------------

    if (model == "theta") {

      lambda <- par[1]
      theta <- par[2]

      yt <- bc(
        y,
        theta
      )

      W <- X

      for (j in seq_len(ncol(W))) {

        if (transform_cols[j]) {

          W[, j] <- bc(
            W[, j],
            lambda
          )
        }
      }

      return(
        list(
          y = yt,
          X = W,
          lambda = lambda,
          theta = theta
        )
      )
    }
  }

  # ============================================================
  # 9. CONCENTRATED LOG-LIKELIHOOD
  # ============================================================

  loglik <- function(par) {

    mod <- construct_model(par)

    fit <- tryCatch(
      lm.fit(
        mod$X,
        mod$y
      ),
      error = function(e)
        NULL
    )

    if (is.null(fit))
      return(-Inf)

    residuals <- fit$residuals

    n <- length(y)

    SSR <- sum(
      residuals^2
    )

    if (!is.finite(SSR) ||
        SSR <= 0) {

      return(-Inf)
    }

    sigma2 <- SSR / n

    ll <-
      -0.5 * n *
      (
        log(2 * pi) +
          1 +
          log(sigma2)
      )

    # ----------------------------------------------------------
    # Jacobian
    # ----------------------------------------------------------

    if (model == "lhsonly") {

      ll <-
        ll +
        (par[1] - 1) *
        sum(log(y))
    }

    if (model == "lambda") {

      ll <-
        ll +
        (par[1] - 1) *
        sum(log(y))
    }

    if (model == "theta") {

      ll <-
        ll +
        (par[2] - 1) *
        sum(log(y))
    }

    ll
  }

  # ============================================================
  # 10. COMPARISON MODEL
  # ============================================================

  comparison_par <-

    if (model == "theta") {

      c(1, 1)

    } else {

      1
    }

  ll_comparison <-
    loglik(
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
      c(1, 0.5)
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
  # 12. OPTIMIZATION
  # ============================================================

  if (trace) {

    cat("\n")
    cat("Fitting comparison model\n\n")

    cat(
      "Iteration 0:  Log likelihood = ",
      formatC(
        ll_comparison,
        digits = 8,
        format = "f"
      ),
      "\n",
      sep = ""
    )

    cat("\n")
    cat("Fitting full model\n\n")
  }

  best <- NULL

  best_value <- Inf

  for (s in seq_len(nrow(starts))) {

    current_iter <- 0

    objective <- function(par) {

      value <- -loglik(par)

      current_iter <<-
        current_iter + 1

      if (trace && s == 1) {

        cat(
          "Iteration ",
          current_iter - 1,
          ":  Log likelihood = ",
          formatC(
            -value,
            digits = 8,
            format = "f"
          ),
          "\n",
          sep = ""
        )
      }

      value
    }

    fit <- tryCatch(
      optim(
        par = starts[s, ],
        fn = objective,
        method = "BFGS",
        hessian = TRUE,
        control = list(
          maxit = 5000,
          reltol = 1e-10
        )
      ),
      error = function(e)
        NULL
    )

    if (!is.null(fit) &&
        is.finite(fit$value)) {

      if (fit$value < best_value) {

        best <- fit

        best_value <-
          fit$value
      }
    }

    # Only show first optimization path
    if (s == 1 && trace)
      cat("\n")
  }

  if (is.null(best)) {

    stop(
      "Maximum likelihood optimization failed."
    )
  }

  parhat <- best$par

  ll_full <-
    -best$value

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

    names(parhat) <-
      c(
        "lambda",
        "theta"
      )
  }

  # ============================================================
  # 14. VARIANCE-COVARIANCE MATRIX
  # ============================================================

  H <- best$hessian

  V <- tryCatch(
    solve(H),
    error = function(e)
      NULL
  )

  if (is.null(V)) {

    H <- tryCatch(
      optimHess(
        parhat,
        function(p)
          -loglik(p)
      ),
      error = function(e)
        NULL
    )

    if (!is.null(H)) {

      V <- tryCatch(
        solve(H),
        error = function(e)
          NULL
      )
    }
  }

  if (!is.null(V)) {

    se_transform <-
      sqrt(
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
  # 15. CONFIDENCE INTERVAL
  # ============================================================

  # IMPORTANT:
  # Stata's level(95), level(99), etc.
  # are percentages, not proportions.

  alpha <-
    1 - level / 100

  zcrit <-
    qnorm(
      1 - alpha / 2
    )

  z_transform <-
    parhat /
    se_transform

  p_transform <-
    2 *
    pnorm(
      abs(z_transform),
      lower.tail = FALSE
    )

  lower <-
    parhat -
    zcrit *
    se_transform

  upper <-
    parhat +
    zcrit *
    se_transform

  transformation_table <-
    data.frame(

      Coefficient =
        parhat,

      `Std. err.` =
        se_transform,

      z =
        z_transform,

      `P>|z|` =
        p_transform,

      Lower =
        lower,

      Upper =
        upper,

      row.names =
        paste0(
          "/",
          names(parhat)
        ),

      check.names =
        FALSE
    )

  # ============================================================
  # 16. SCALE-VARIANT PARAMETERS
  # ============================================================

  mod_full <-
    construct_model(
      parhat
    )

  fit_full <-
    lm.fit(
      mod_full$X,
      mod_full$y
    )

  beta <-
    fit_full$coefficients

  sigma <-
    sqrt(
      sum(
        fit_full$residuals^2
      ) / length(y)
    )

  names(beta) <-
    colnames(
      mod_full$X
    )

  # ============================================================
  # 17. LR TESTS
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

  LR_tests <- list()

  for (nm in names(restricted_values)) {

    ll0 <-
      loglik(
        restricted_values[[nm]]
      )

    LR <-
      max(
        0,
        2 *
          (
            ll_full -
              ll0
          )
      )

    df_lr <-
      length(parhat)

    p_lr <-
      pchisq(
        LR,
        df = df_lr,
        lower.tail = FALSE
      )

    LR_tests[[nm]] <-
      c(
        Restricted =
          ll0,

        `LR statistic` =
          LR,

        df =
          df_lr,

        `Prob > chi2` =
          p_lr
      )
  }

  LR_table <-
    do.call(
      rbind,
      LR_tests
    )

  test_parameter <-

    if (model == "lhsonly" ||
        model == "theta") {

      "theta"

    } else {

      "lambda"
    }

  rownames(LR_table) <-
    paste0(
      test_parameter,
      " = ",
      names(restricted_values)
    )

  # ============================================================
  # 18. MODEL LR TEST
  # ============================================================

  LR_model <-
    max(
      0,
      2 *
        (
          ll_full -
            ll_comparison
        )
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
  # 19. FORMATTING
  # ============================================================

  fmt <- function(x) {

    ifelse(
      is.na(x),
      "NA",
      formatC(
        x,
        digits = digits,
        format = "f"
      )
    )
  }

  fmt_p <- function(x) {

    ifelse(
      x < 0.001,
      "0.000",
      fmt(x)
    )
  }

  # ============================================================
  # 20. HEADER
  # ============================================================

  cat(
    "Number of obs   = ",
    nrow(mf),
    "\n",
    sep = ""
  )

  cat(
    "LR chi2(",
    df_model,
    ")      = ",
    formatC(
      LR_model,
      digits = 2,
      format = "f"
    ),
    "\n",
    sep = ""
  )

  cat(
    "Log likelihood = ",
    fmt(ll_full),
    "\n",
    sep = ""
  )

  cat(
    "Prob > chi2     = ",
    fmt_p(p_model),
    "\n\n",
    sep = ""
  )

  # ============================================================
  # 21. TRANSFORMATION TABLE
  # ============================================================

  cat(
    "------------------------------------------------------------------------------\n"
  )

  dep_name <-
    deparse(
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

  for (i in seq_len(
    nrow(
      transformation_table
    )
  )) {

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

        fmt_p(
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
  # 22. SCALE-VARIANT PARAMETERS
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
  # 23. LR TESTS
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

  for (i in seq_len(
    nrow(LR_table)
  )) {

    cat(
      sprintf(
        "%12s %16s %12s %15s\n",

        rownames(
          LR_table
        )[i],

        fmt(
          LR_table[
            i,
            "Restricted"
          ]
        ),

        fmt(
          LR_table[
            i,
            "LR statistic"
          ]
        ),

        fmt_p(
          LR_table[
            i,
            "Prob > chi2"
          ]
        )
      )
    )
  }

  cat(
    "---------------------------------------------------------\n"
  )

  # ============================================================
  # 24. RETURN OBJECT
  # ============================================================

  result <- list(

    call =
      match.call(),

    model =
      model,

    n =
      nrow(mf),

    logLik =
      ll_full,

    comparison_logLik =
      ll_comparison,

    LR =
      LR_model,

    df =
      df_model,

    p_value =
      p_model,

    transformation =
      transformation_table,

    coefficients =
      beta,

    sigma =
      sigma,

    LR_tests =
      LR_table,

    vcov =
      V,

    hessian =
      H,

    convergence =
      best$convergence,

    iterations =
      best$counts["function"]
  )

  class(result) <-
    "boxcox_r"

  invisible(result)
}
