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
#' @param noconstant Suppress the constant term.
#' @param trace Display optimization iterations.
#'
#' @return A Box-Cox regression result.
#' @export

boxcox <- function(formula,
                   data,
                   model = c(
                     "lhsonly",
                     "rhsonly",
                     "lambda",
                     "theta"
                   ),
                   notrans = NULL,
                   level = 95,
                   digits = 4,
                   noconstant = FALSE,
                   trace = TRUE) {

  # ============================================================
  # 1. VALIDATIONS
  # ============================================================

  model <- match.arg(model)

  if (!is.numeric(level) ||
      length(level) != 1 ||
      level <= 0 ||
      level >= 100) {

    stop(
      "'level' must be between 0 and 100."
    )
  }

  if (!is.numeric(digits) ||
      length(digits) != 1 ||
      digits < 0 ||
      digits != floor(digits)) {

    stop(
      "'digits' must be a non-negative integer."
    )
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

    stop(
      "The dependent variable must be numeric."
    )
  }

  if (any(!is.finite(y))) {

    stop(
      "The dependent variable contains non-finite values."
    )
  }

  # ============================================================
  # 3. MODEL TERMS
  # ============================================================

  terms_obj <- terms(formula)

  rhs_vars <- all.vars(
    delete.response(
      terms_obj
    )
  )

  # Variables appearing in the RHS
  rhs_terms <- attr(
    terms_obj,
    "term.labels"
  )

  # ============================================================
  # 4. NOTRANS
  # ============================================================

  if (is.null(notrans)) {

    notrans <- character(0)

  }

  if (!is.character(notrans)) {

    stop(
      "'notrans' must contain variable names."
    )
  }

  unknown_notrans <- setdiff(
    notrans,
    rhs_vars
  )

  if (length(unknown_notrans) > 0) {

    stop(
      paste0(
        "Variables in 'notrans' not found in the RHS: ",
        paste(
          unknown_notrans,
          collapse = ", "
        )
      )
    )
  }

  # ============================================================
  # 5. VARIABLES TO TRANSFORM
  # ============================================================

  transform_x <-
    model %in% c(
      "rhsonly",
      "lambda",
      "theta"
    )

  transformed_vars <- if (transform_x) {

    setdiff(
      rhs_vars,
      notrans
    )

  } else {

    character(0)

  }

  # ============================================================
  # 6. CHECK POSITIVITY
  # ============================================================

  # ------------------------------------------------------------
  # LHS
  # ------------------------------------------------------------

  if (model %in% c(
    "lhsonly",
    "lambda",
    "theta"
  )) {

    if (any(y <= 0)) {

      stop(
        "The dependent variable must be strictly positive ",
        "for this Box-Cox model."
      )
    }
  }

  # ------------------------------------------------------------
  # RHS
  # ------------------------------------------------------------

  if (transform_x &&
      length(transformed_vars) > 0) {

    for (v in transformed_vars) {

      if (!v %in% names(mf))
        next

      xv <- mf[[v]]

      if (!is.numeric(xv)) {

        stop(
          paste0(
            "Variable '",
            v,
            "' must be numeric because it is transformed."
          )
        )
      }

      if (any(!is.finite(xv))) {

        stop(
          paste0(
            "Variable '",
            v,
            "' contains non-finite values."
          )
        )
      }

      if (any(xv <= 0)) {

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
  # 7. DESIGN MATRIX
  # ============================================================

  X <- model.matrix(
    formula,
    mf
  )

  if (noconstant &&
      "(Intercept)" %in% colnames(X)) {

    X <- X[
      ,
      colnames(X) != "(Intercept)",
      drop = FALSE
    ]
  }

  # ============================================================
  # 8. IDENTIFY RHS COLUMNS
  # ============================================================

  # Map each column of model.matrix() to its formula term.
  assign_vector <- attr(
    X,
    "assign"
  )

  # Remove intercept from the transformation candidates
  candidate_cols <- seq_len(
    ncol(X)
  )

  if ("(Intercept)" %in% colnames(X)) {

    candidate_cols <- candidate_cols[
      colnames(X) != "(Intercept)"
    ]
  }

  transform_cols <- rep(
    FALSE,
    ncol(X)
  )

  # ------------------------------------------------------------
  # Simple-variable formulas
  # ------------------------------------------------------------

  for (v in transformed_vars) {

    # Find exact RHS term
    term_index <- which(
      rhs_terms == v
    )

    if (length(term_index) == 1) {

      transform_cols[
        assign_vector == term_index
      ] <- TRUE
    }
  }

  # ============================================================
  # 9. BOX-COX TRANSFORMATION
  # ============================================================

  bc <- function(x, lambda) {

    lx <- log(x)

    if (abs(lambda) < 1e-8) {

      return(lx)

    }

    expm1(
      lambda * lx
    ) / lambda
  }

  # ============================================================
  # 10. RHS JACOBIAN
  # ============================================================

  rhs_log_jacobian <- function(lambda) {

    if (!transform_x ||
        length(transformed_vars) == 0) {

      return(0)
    }

    total <- 0

    for (v in transformed_vars) {

      xv <- mf[[v]]

      total <-
        total +
        (lambda - 1) *
        sum(log(xv))
    }

    total
  }

  # ============================================================
  # 11. CONSTRUCT MODEL
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

          # Original variable corresponding to this column
          term_index <-
            assign_vector[j]

          term_name <-
            rhs_terms[term_index]

          if (term_name %in% transformed_vars) {

            W[, j] <-
              bc(
                W[, j],
                lambda
              )
          }
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

          term_index <-
            assign_vector[j]

          term_name <-
            rhs_terms[term_index]

          if (term_name %in% transformed_vars) {

            W[, j] <-
              bc(
                W[, j],
                lambda
              )
          }
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

          term_index <-
            assign_vector[j]

          term_name <-
            rhs_terms[term_index]

          if (term_name %in% transformed_vars) {

            W[, j] <-
              bc(
                W[, j],
                lambda
              )
          }
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
  # 12. CONCENTRATED LOG-LIKELIHOOD
  # ============================================================

  loglik <- function(par) {

    mod <- construct_model(
      par
    )

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

    sigma2 <-
      SSR / n

    ll <-
      -0.5 *
      n *
      (
        log(2 * pi) +
          1 +
          log(sigma2)
      )

    # ----------------------------------------------------------
    # LHS Jacobian
    # ----------------------------------------------------------

    if (model == "lhsonly") {

      theta <- par[1]

      ll <-
        ll +
        (theta - 1) *
        sum(log(y))
    }

    if (model == "lambda") {

      lambda <- par[1]

      ll <-
        ll +
        (lambda - 1) *
        sum(log(y))
    }

    if (model == "theta") {

      theta <- par[2]

      ll <-
        ll +
        (theta - 1) *
        sum(log(y))
    }

    # ----------------------------------------------------------
    # RHS Jacobian
    # ----------------------------------------------------------

    if (model == "rhsonly") {

      lambda <- par[1]

      ll <-
        ll +
        rhs_log_jacobian(
          lambda
        )
    }

    if (model == "lambda") {

      lambda <- par[1]

      ll <-
        ll +
        rhs_log_jacobian(
          lambda
        )
    }

    if (model == "theta") {

      lambda <- par[1]

      ll <-
        ll +
        rhs_log_jacobian(
          lambda
        )
    }

    ll
  }

  # ============================================================
  # 13. COMPARISON MODEL
  # ============================================================

  comparison_par <-

    if (model == "theta") {

      c(
        1,
        1
      )

    } else {

      1
    }

  ll_comparison <-
    loglik(
      comparison_par
    )

  # ============================================================
  # 14. STARTING VALUES
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
  # 15. OPTIMIZATION
  # ============================================================

  if (trace) {

    cat(
      "\nFitting comparison model\n\n"
    )

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

    cat(
      "\nFitting full model\n\n"
    )
  }

  best <- NULL

  best_value <- Inf

  for (s in seq_len(
    nrow(starts)
  )) {

    current_iter <- 0

    objective <- function(par) {

      value <- -loglik(
        par
      )

      current_iter <<-
        current_iter + 1

      if (trace &&
          s == 1) {

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

    if (s == 1 &&
        trace) {

      cat("\n")
    }
  }

  if (is.null(best)) {

    stop(
      "Maximum likelihood optimization failed."
    )
  }

  # ============================================================
  # 16. ESTIMATES
  # ============================================================

  parhat <- best$par

  ll_full <-
    -best$value

  if (model == "lhsonly") {

    names(parhat) <-
      "theta"

  } else if (model == "rhsonly") {

    names(parhat) <-
      "lambda"

  } else if (model == "lambda") {

    names(parhat) <-
      "lambda"

  } else {

    names(parhat) <-
      c(
        "lambda",
        "theta"
      )
  }

  # ============================================================
  # 17. VARIANCE-COVARIANCE MATRIX
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
  # 18. CONFIDENCE INTERVAL
  # ============================================================

  alpha <-
    1 -
    level / 100

  zcrit <-
    qnorm(
      1 -
        alpha / 2
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

      check.names = FALSE
    )

  # ============================================================
  # 19. SCALE-VARIANT PARAMETERS
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

  names(beta) <-
    colnames(
      mod_full$X
    )

  sigma <-
    sqrt(
      sum(
        fit_full$residuals^2
      ) /
        length(y)
    )

  # ============================================================
  # 20. CLASSIFY SCALE-VARIANT PARAMETERS
  # ============================================================

  beta_notrans <- character(0)

  beta_trans <- character(0)

  if ("(Intercept)" %in% names(beta)) {

    beta_notrans <-
      c(
        beta_notrans,
        "(Intercept)"
      )
  }

  # Explicit notrans variables
  if (length(notrans) > 0) {

    for (v in notrans) {

      idx <- which(
        names(beta) == v
      )

      if (length(idx) == 1) {

        beta_notrans <-
          c(
            beta_notrans,
            v
          )
      }
    }
  }

  # Transformed variables
  if (transform_x &&
      length(transformed_vars) > 0) {

    for (v in transformed_vars) {

      idx <- which(
        names(beta) == v
      )

      if (length(idx) == 1) {

        beta_trans <-
          c(
            beta_trans,
            v
          )
      }
    }
  }

  # ============================================================
  # 21. LR TESTS OF TRANSFORMATION PARAMETERS
  # ============================================================

  restricted_values <-

    if (model == "theta") {

      list(

        "-1" =
          c(-1, -1),

        "0" =
          c(0, 0),

        "1" =
          c(1, 1)

      )

    } else {

      list(

        "-1" =
          -1,

        "0" =
          0,

        "1" =
          1

      )
    }

  LR_tests <- list()

  for (nm in names(
    restricted_values
  )) {

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

    if (model %in% c(
      "lhsonly",
      "theta"
    )) {

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
  # 22. MODEL LR TEST
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

  # Stata's displayed LR chi2 degrees of freedom
  # correspond to the number of RHS model terms.
  df_model <-
    length(
      rhs_terms
    )

  p_model <-
    pchisq(
      LR_model,
      df = df_model,
      lower.tail = FALSE
    )

  # ============================================================
  # 23. FORMATTING
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
  # 24. HEADER
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
  # 25. TRANSFORMATION TABLE
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
  # 26. SCALE-VARIANT PARAMETERS
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

  # ------------------------------------------------------------
  # NOTRANS
  # ------------------------------------------------------------

  cat(
    "Notrans      |\n"
  )

  if (length(beta_notrans) > 0) {

    for (nm in beta_notrans) {

      display_name <-
        ifelse(
          nm == "(Intercept)",
          "_cons",
          nm
        )

      cat(
        sprintf(
          "%12s | %12s\n",
          display_name,
          fmt(
            beta[nm]
          )
        )
      )
    }
  }

  # ------------------------------------------------------------
  # TRANS
  # ------------------------------------------------------------

  if (length(beta_trans) > 0) {

    cat(
      "-------------+--------------\n"
    )

    cat(
      "Trans        |\n"
    )

    for (nm in beta_trans) {

      cat(
        sprintf(
          "%12s | %12s\n",
          nm,
          fmt(
            beta[nm]
          )
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
  # 27. LR TEST TABLE
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
  # 28. RETURN OBJECT
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
      best$counts[
        "function"
      ],

    transformed_vars =
      transformed_vars,

    notrans =
      notrans
  )

  class(result) <-
    "boxcox_r"

  invisible(result)
}
