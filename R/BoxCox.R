#' Stata-like Box-Cox Regression
#'
#' @description
#' Estimates Box-Cox regression models using Maximum Likelihood, replicating the exact 
#' behavior, scaling, and output format of Stata's `boxcox` command. It supports 
#' simultaneous transformation parameters (theta and lambda) and allows specifying 
#' variables to remain untransformed.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class).
#' @param data a data frame, list or environment containing the variables in the model.
#' @param notrans a character vector specifying the names of independent variables that 
#' should NOT be transformed (e.g., dummy variables).
#' @param model character string specifying the model type: \code{"theta"} (default, 
#' different parameters for LHS and RHS), \code{"lambda"} (same parameter), 
#' \code{"lhsonly"}, or \code{"rhsonly"}.
#' @param init_vals numeric vector of initial values for the optimizer.
#' @param level numeric value between 0 and 1 specifying the confidence level. Default is 0.95.
#' @param digits integer indicating the number of decimal places to print. Default is 6.
#' @param stata_df_compat logical. If \code{TRUE} (default), replicates Stata's historical 
#' behavior of not subtracting omitted collinear variables from the global Chi2 degrees of 
#' freedom, and forcing 1 df for the lower restricted tests. If \code{FALSE}, calculates 
#' statistically correct degrees of freedom.
#' @param trace logical. If \code{TRUE} (default), prints the iteration log for the 
#' likelihood maximization process.
#'
#' @return A list containing the estimated parameters, log-likelihood, coefficients, and standard errors.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'df' is your dataset and 'dummy1' is a binary variable
#' result <- box.cox(y ~ x1 + x2 + dummy1, data = df, 
#'                   notrans = "dummy1", model = "theta")
#' }
box.cox <- function(formula, data, notrans = NULL, model = "theta", 
                    init_vals = NULL, level = 0.95, digits = 6,
                    stata_df_compat = TRUE, trace = TRUE) {
  
  # -------------------------------------------------------------
  # 1. INPUT VALIDATION & SETUP
  # -------------------------------------------------------------
  if (!inherits(formula, "formula")) stop("'formula' must be of class formula.")
  if (missing(data) || !is.data.frame(data)) stop("'data' must be a data.frame.")
  
  valid_models <- c("theta", "lambda", "lhsonly", "rhsonly")
  if (!(model %in% valid_models)) stop("Model must be one of: ", paste(valid_models, collapse = ", "))
  
  if (level > 1) level <- level / 100 
  if (level <= 0 || level >= 1) stop("'level' must be between 0 and 1 (e.g., 0.95).")
  z_crit <- abs(qnorm((1 - level) / 2))
  
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  y_name <- names(mf)[1] 
  
  X_mat <- model.matrix(formula, mf)
  if ("(Intercept)" %in% colnames(X_mat)) {
    X_mat <- X_mat[, colnames(X_mat) != "(Intercept)", drop = FALSE]
  }
  
  if (any(y <= 0, na.rm = TRUE)) stop("The dependent variable (", y_name, ") must be strictly positive.")
  
  if (!is.null(notrans)) {
    missing_vars <- setdiff(notrans, colnames(X_mat))
    if (length(missing_vars) > 0) {
      stop("Variables specified in 'notrans' not found in formula: ", paste(missing_vars, collapse = ", "))
    }
    
    X_notrans <- X_mat[, notrans, drop = FALSE]
    trans_cols <- setdiff(colnames(X_mat), notrans)
    if (length(trans_cols) == 0 && model != "lhsonly") {
      stop("At least one independent variable must remain for transformation.")
    }
    X_trans <- X_mat[, trans_cols, drop = FALSE]
  } else {
    X_trans <- X_mat
    X_notrans <- NULL
  }
  
  if (model != "lhsonly") {
    if (any(X_trans <= 0, na.rm = TRUE)) stop("Independent variables to transform must be strictly positive.")
  }
  
  var_names_trans <- colnames(X_trans)
  var_names_notrans <- if (!is.null(X_notrans)) colnames(X_notrans) else NULL
  N <- length(y)
  
  bc_std <- function(v, lambda) {
    if (abs(lambda) < 1e-5) log(v) else (v^lambda - 1) / lambda
  }
  
  calc_ll <- function(t_val, l_val) {
    yt <- bc_std(y, t_val)
    Xt <- if (model == "lhsonly") X_trans else apply(X_trans, 2, function(col) bc_std(col, l_val))
    
    X_all <- if (is.null(X_notrans)) cbind(1, Xt) else cbind(1, Xt, X_notrans)
    
    if (any(!is.finite(yt)) || any(!is.finite(X_all))) return(-1e10)
    res <- tryCatch(qr.resid(qr(X_all), yt), error = function(e) NULL)
    if (is.null(res)) return(-1e10)
    
    sigma2 <- sum(res^2) / N
    if (is.na(sigma2) || sigma2 <= 0) return(-1e10)
    
    return(- (N/2)*log(2*pi) - (N/2)*log(sigma2) - (N/2) + (t_val - 1)*sum(log(y)))
  }
  
  # -------------------------------------------------------------
  # 2. OPTIMIZATION WITH ITERATION TRACKER
  # -------------------------------------------------------------
  if (is.null(init_vals)) init_vals <- if(model == "theta") c(1, 1) else 1
  
  get_params <- function(params) {
    if (model == "lhsonly") return(c(params[1], 1))
    if (model == "rhsonly") return(c(1, params[1]))
    if (model == "lambda")  return(c(params[1], params[1]))
    return(c(params[1], params[2]))
  }
  
  run_optim <- function(start_vals, title) {
    tracker <- new.env()
    tracker$best_ll <- -Inf
    tracker$iter <- 0
    
    if (trace) cat(paste0("\n", title, "\n\n"))
    
    obj_fun <- function(params) {
      p <- get_params(params)
      current_ll <- calc_ll(p[1], p[2])
      
      if (trace && current_ll > (tracker$best_ll + 1e-5)) {
        cat(sprintf("Iteration %-2d:  Log likelihood = %10.6f\n", tracker$iter, current_ll))
        tracker$best_ll <- current_ll
        tracker$iter <- tracker$iter + 1
      }
      return(-current_ll)
    }
    
    opt <- optim(par = start_vals, fn = obj_fun, method = "L-BFGS-B", 
                 lower = -5, upper = 5, hessian = TRUE,
                 control = list(maxit = 3000, factr = 1e7))
    return(opt)
  }
  
  opt_full <- run_optim(init_vals, "Fitting full model")
  p_final <- get_params(opt_full$par)
  t_final <- p_final[1]
  l_final <- p_final[2]
  ll_full <- -opt_full$value
  
  # -------------------------------------------------------------
  # 3. POST-ESTIMATION COMPUTATIONS
  # -------------------------------------------------------------
  vcov_mat <- tryCatch(solve(opt_full$hessian), error = function(e) matrix(NA, nrow=length(opt_full$par), ncol=length(opt_full$par)))
  se_params <- sqrt(abs(diag(vcov_mat)))
  z_vals <- opt_full$par / se_params
  p_vals <- 2 * (1 - pnorm(abs(z_vals)))
  ci_low <- opt_full$par - z_crit * se_params
  ci_high <- opt_full$par + z_crit * se_params
  
  yt_final <- bc_std(y, t_final)
  Xt_final <- if (model == "lhsonly") X_trans else apply(X_trans, 2, function(col) bc_std(col, l_final))
  
  X_final_mat <- cbind(Intercept = 1, Xt_final)
  colnames(X_final_mat) <- c("_cons", var_names_trans)
  if (!is.null(X_notrans)) X_final_mat <- cbind(X_final_mat, X_notrans)
  
  fit_final <- lm.fit(X_final_mat, yt_final)
  betas <- fit_final$coefficients
  omitidas <- names(betas)[is.na(betas)]
  sigma_val <- sqrt(sum(fit_final$residuals^2) / N)
  
  k_inicial <- length(var_names_trans) + length(var_names_notrans)
  k_eff <- length(betas) - 1 - length(omitidas)
  df_global_base <- if (model %in% c("theta", "rhsonly")) 1 else 0
  
  if (stata_df_compat) {
    df_global <- k_inicial + df_global_base
    df_test <- 1 
  } else {
    df_global <- k_eff + df_global_base 
    df_test <- if (model == "theta") 2 else 1 
  }
  
  if (model == "rhsonly") {
    res_comp <- y - mean(y)
    ll_comp <- - (N/2)*log(2*pi) - (N/2)*log(sum(res_comp^2)/N) - (N/2)
  } else {
    obj_comp <- function(th) {
      yt_c <- bc_std(y, th)
      s2_c <- sum((yt_c - mean(yt_c))^2) / N
      if (is.na(s2_c) || s2_c <= 0) return(1e10)
      return(- (- (N/2)*log(2*pi) - (N/2)*log(s2_c) - (N/2) + (th - 1)*sum(log(y))))
    }
    if (trace) cat(paste0("\nFitting comparison model\n\n"))
    tracker_comp <- new.env(); tracker_comp$best_ll <- -Inf; tracker_comp$iter <- 0
    obj_comp_tracked <- function(th) {
      val <- -obj_comp(th)
      if (trace && val > (tracker_comp$best_ll + 1e-5)) {
        cat(sprintf("Iteration %-2d:  Log likelihood = %10.6f\n", tracker_comp$iter, val))
        tracker_comp$best_ll <- val
        tracker_comp$iter <- tracker_comp$iter + 1
      }
      return(-val)
    }
    opt_comp <- optim(par = 1, fn = obj_comp_tracked, method = "L-BFGS-B", lower = -5, upper = 5)
    ll_comp <- -opt_comp$value
  }
  
  lr_chi_main <- 2 * max(0, ll_full - ll_comp)
  p_chi_main  <- pchisq(lr_chi_main, df = df_global, lower.tail = FALSE)
  
  get_test_ll <- function(val) {
    if (model == "lhsonly") return(calc_ll(val, 1))
    if (model == "rhsonly") return(calc_ll(1, val))
    return(calc_ll(val, val))
  }
  
  test_label <- switch(model, "lhsonly" = "theta = ", "rhsonly" = "lambda = ", "theta=lambda = ")
  
  # -------------------------------------------------------------
  # 4. STATA-STYLE CONSOLE OUTPUT
  # -------------------------------------------------------------
  cat("\n")
  if (length(omitidas) > 0) {
    cat(sprintf("note: %s omitted because of collinearity.\n\n", paste(omitidas, collapse = ", ")))
  }
  
  fmt_ll <- paste0("%-15.", digits, "f")
  
  cat(sprintf("%50s = %10d\n", "Number of obs", N))
  cat(sprintf("%50s = %10.2f\n", paste0("LR chi2(", df_global, ")"), max(0, lr_chi_main)))
  cat(sprintf(paste0("Log likelihood = ", fmt_ll, "%24s = %10.3f\n"), ll_full, "Prob > chi2", p_chi_main))
  
  cat(strrep("-", 78), "\n")
  y_str <- substr(y_name, 1, 13)
  
  conf_str <- paste0("[", round(level * 100, 1), "% conf. interval]")
  cat(sprintf("%13s | %10s  %10s %5s %8s   %s\n", y_str, "Coefficient", "Std. err.", "z", "P>|z|", conf_str))
  cat(strrep("-", 13), "+", strrep("-", 62), "\n", sep="")
  
  print_row <- function(name, coef, se, z, p, ci_l, ci_h) {
    fmt_row <- paste0("%13s | %10.", digits, "f  %10.2e %5.2f %8.3f   %11.", digits, "f  %11.", digits, "f\n")
    cat(sprintf(fmt_row, name, coef, se, z, p, ci_l, ci_h))
  }
  
  if (model %in% c("lambda", "rhsonly")) {
    print_row("/lambda", opt_full$par[1], se_params[1], z_vals[1], p_vals[1], ci_low[1], ci_high[1])
  } else if (model == "lhsonly") {
    print_row("/theta", opt_full$par[1], se_params[1], z_vals[1], p_vals[1], ci_low[1], ci_high[1])
  } else if (model == "theta") {
    print_row("/lambda", l_final, se_params[2], z_vals[2], p_vals[2], ci_low[2], ci_high[2])
    print_row("/theta", t_final, se_params[1], z_vals[1], p_vals[1], ci_low[1], ci_high[1])
  }
  cat(strrep("-", 78), "\n\n")
  
  cat("Estimates of scale-variant parameters\n")
  cat(strrep("-", 28), "\n")
  cat(sprintf("%13s| %12s\n", " ", "Coefficient"))
  cat(strrep("-", 13), "+", strrep("-", 14), "\n", sep="")
  
  fmt_scale <- paste0("%13s| %12.", digits, "g\n")
  
  cat(sprintf("%-13s|\n", "Notrans"))
  cat(sprintf(fmt_scale, "_cons", betas["_cons"]))
  
  if (model == "lhsonly") {
    for (v in var_names_trans) if (!(v %in% omitidas)) cat(sprintf(fmt_scale, v, betas[v]))
  }
  if (!is.null(X_notrans)) {
    for (v in var_names_notrans) if (!(v %in% omitidas)) cat(sprintf(fmt_scale, v, betas[v]))
  }
  cat(strrep("-", 13), "+", strrep("-", 14), "\n", sep="")
  
  if (model != "lhsonly") {
    cat(sprintf("%-13s|\n", "Trans"))
    for (v in var_names_trans) if (!(v %in% omitidas)) cat(sprintf(fmt_scale, v, betas[v]))
    cat(strrep("-", 13), "+", strrep("-", 14), "\n", sep="")
  }
  cat(sprintf(fmt_scale, "/sigma", sigma_val))
  cat(strrep("-", 28), "\n\n")
  
  cat(strrep("-", 63), "\n")
  cat(sprintf("   %-15s    %-15s\n", "Test", "Restricted"))
  cat(sprintf("   %-15s    %-15s %10s %15s\n", "H0:", "log likelihood", "chi2", "Prob > chi2"))
  cat(strrep("-", 63), "\n")
  
  print_test <- function(val_str, ll_rest) {
    chi2 <- 2 * max(0, ll_full - ll_rest)
    p_val <- pchisq(chi2, df = df_test, lower.tail = FALSE)
    fmt_test <- paste0("%-20s %10.", digits, "f       %10.2f %15.3f\n")
    cat(sprintf(fmt_test, paste0(test_label, val_str), ll_rest, chi2, p_val))
  }
  
  print_test("-1", get_test_ll(-1))
  print_test(" 0", get_test_ll(0))
  print_test(" 1", get_test_ll(1))
  cat(strrep("-", 63), "\n\n")
  
  res_list <- list(
    theta = t_final, lambda = l_final, loglik = ll_full, 
    coefficients = betas, sigma = sigma_val, 
    vcov = vcov_mat, p_values = p_vals
  )
  class(res_list) <- "box.cox"
  return(invisible(res_list))
}
