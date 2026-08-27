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

boxcox <- function(formula, data, model = "lhsonly", notrans = NULL,
                   level = 0.95) {
  
  model <- match.arg(model, c("lhsonly", "rhsonly", "lambda", "theta"))
  
  mf <- model.frame(formula, data = data, na.action = na.omit)
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  
  if (any(y <= 0))
    stop("La variable dependiente debe ser estrictamente positiva.")
  
  xnames <- colnames(X)
  xnames <- xnames[xnames != "(Intercept)"]
  
  if (!is.null(notrans)) {
    trans <- !(xnames %in% notrans)
  } else {
    trans <- rep(TRUE, length(xnames))
  }
  
  Xtrans <- X[, c(FALSE, trans), drop = FALSE]
  Xnotrans <- X[, c(FALSE, !trans), drop = FALSE]
  
  if (ncol(Xtrans) > 0 && any(Xtrans <= 0))
    stop("Las variables independientes transformadas deben ser estrictamente positivas.")
  
  bc <- function(x, lambda) {
    if (abs(lambda) < 1e-8)
      log(x)
    else
      (x^lambda - 1) / lambda
  }
  
  loglik <- function(par) {
    
    if (model == "lhsonly") {
      theta <- par[1]
      yt <- bc(y, theta)
      W <- X
      
    } else if (model == "rhsonly") {
      lambda <- par[1]
      Xt <- if (ncol(Xtrans) > 0)
        apply(Xtrans, 2, bc, lambda = lambda)
      else NULL
      
      W <- cbind(1, Xt, Xnotrans)
      yt <- y
      
    } else if (model == "lambda") {
      lambda <- par[1]
      
      yt <- bc(y, lambda)
      
      Xt <- if (ncol(Xtrans) > 0)
        apply(Xtrans, 2, bc, lambda = lambda)
      else NULL
      
      W <- cbind(1, Xt, Xnotrans)
      
    } else {
      lambda <- par[1]
      theta <- par[2]
      
      yt <- bc(y, theta)
      
      Xt <- if (ncol(Xtrans) > 0)
        apply(Xtrans, 2, bc, lambda = lambda)
      else NULL
      
      W <- cbind(1, Xt, Xnotrans)
    }
    
    fit <- lm.fit(W, yt)
    e <- fit$residuals
    n <- length(y)
    sigma2 <- sum(e^2) / n
    
    if (!is.finite(sigma2) || sigma2 <= 0)
      return(-Inf)
    
    jac <- 0
    
    if (model %in% c("lhsonly", "lambda", "theta")) {
      theta_y <- if (model == "lambda") par[1] else
        if (model == "theta") par[2] else par[1]
      
      jac <- (theta_y - 1) * sum(log(y))
    }
    
    (-n / 2) * (log(2 * pi) + 1 + log(sigma2)) + jac
  }
  
  start <- switch(
    model,
    lhsonly = 1,
    rhsonly = 1,
    lambda = 1,
    theta = c(1, 1)
  )
  
  opt <- optim(
    start,
    function(p) -loglik(p),
    method = "BFGS",
    hessian = TRUE
  )
  
  est <- opt$par
  ll_full <- -opt$value
  
  if (model == "lhsonly")
    names(est) <- "theta"
  
  if (model == "rhsonly")
    names(est) <- "lambda"
  
  if (model == "lambda")
    names(est) <- "lambda"
  
  if (model == "theta")
    names(est) <- c("lambda", "theta")
  
  # Error estándar
  V <- tryCatch(
    solve(opt$hessian),
    error = function(e) matrix(NA, length(est), length(est))
  )
  
  se <- sqrt(diag(V))
  
  # Wald z y p-value
  z <- est / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  zcrit <- qnorm(1 - (1 - level) / 2)
  
  lower <- est - zcrit * se
  upper <- est + zcrit * se
  
  results <- data.frame(
    Estimate = est,
    Std_Error = se,
    z = z,
    p_value = p,
    Lower = lower,
    Upper = upper
  )
  
  rownames(results) <- names(est)
  
  # LR tests
  lr_test <- function(restricted) {
    
    ll_r <- loglik(restricted)
    LR <- 2 * (ll_full - ll_r)
    
    df <- length(est)
    pval <- pchisq(LR, df = df, lower.tail = FALSE)
    
    c(LR = LR, df = df, p_value = pval)
  }
  
  if (model == "theta") {
    test_1 <- lr_test(c(1, 1))
    test_0 <- lr_test(c(0, 0))
    test_m1 <- lr_test(c(-1, -1))
  } else {
    test_1 <- lr_test(1)
    test_0 <- lr_test(0)
    test_m1 <- lr_test(-1)
  }
  
  LRtests <- rbind(
    "lambda/theta = 1" = test_1,
    "lambda/theta = 0" = test_0,
    "lambda/theta = -1" = test_m1
  )
  
  # Salida
  cat("\nBox-Cox regression model\n")
  cat("Model:", model, "\n\n")
  
  cat("Transformation parameters:\n\n")
  
  print(
    round(results, 4),
    digits = 4
  )
  
  cat("\nLikelihood-ratio tests:\n\n")
  
  print(
    round(LRtests, 4),
    digits = 4
  )
  
  cat("\nLog likelihood =", round(ll_full, 4), "\n")
  
  invisible(
    list(
      model = model,
      coefficients = results,
      LR_tests = LRtests,
      logLik = ll_full,
      estimates = est,
      vcov = V
    )
  )
}
