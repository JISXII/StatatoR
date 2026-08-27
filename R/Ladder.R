#' Ladder of Powers Normality Tests (All Methods: Stata, Shapiro, Jarque-Bera)
#'
#' @description
#' Searches the ladder of powers (Tukey's transformations) to convert a variable into 
#' a normally distributed one. It outputs three distinct categorized tables in the console: 
#' Stata's Skewness-Kurtosis test (sktest), Shapiro-Wilk test, and Jarque-Bera test.
#'
#' @param x A numeric vector to be transformed.
#' @param na.rm Logical. Should missing values be removed? Default is TRUE.
#'
#' @return A list containing N and the test statistics for each transformation, returned invisibly.
#' @export
ladder <- function(x, na.rm = TRUE) {
  
  if (!is.numeric(x)) stop("'x' must be a numeric vector.")
  if (na.rm) x <- na.omit(x)
  
  if (any(x <= 0)) {
    message("note: variable contains zero or negative values. Log and root transformations will produce NAs.")
  }
  
  n_total <- length(x)
  if (n_total < 8) stop("number of observations must be at least 8 for normality checks.")
  
  # -------------------------------------------------------------
  # 1. MÉTODOS ESTADÍSTICOS
  # -------------------------------------------------------------
  
  # A. Método Stata (sktest: D'Agostino-Pearson con ajuste de Royston)[cite: 2]
  stata_sktest <- function(v) {
    n <- length(v)
    m1 <- mean(v)
    m2 <- sum((v - m1)^2) / n
    if (m2 == 0) return(c(chi2 = NA, p = NA))
    m3 <- sum((v - m1)^3) / n
    m4 <- sum((v - m1)^4) / n
    
    g1 <- m3 / (m2^(3/2))
    y_skew <- g1 * sqrt(((n + 1) * (n + 3)) / (6 * (n - 2)))
    beta2 <- (3 * (n^2 + 27*n - 70) * (n + 1) * (n + 3)) / ((n - 2) * (n + 5) * (n + 7) * (n + 9))
    w2 <- sqrt(2 * beta2 - 2) - 1
    delta <- 1 / sqrt(log(sqrt(w2)))
    alpha <- sqrt(2 / (w2 - 1))
    z1 <- delta * log(y_skew / alpha + sqrt((y_skew / alpha)^2 + 1))
    
    g2 <- (m4 / (m2^2)) - 3
    mean_g2 <- (6 * (n - 2)) / ((n + 1) * (n + 3))
    var_g2 <- (24 * n * (n - 2) * (n - 3)) / (((n + 1)^2) * (n + 3) * (n + 5))
    x_kurt <- (g2 - mean_g2) / sqrt(var_g2)
    
    s_kurt <- (6 * (n^2 - 5*n + 2)) / ((n + 7) * (n + 9)) * sqrt((6 * (n + 3) * (n + 5)) / (n * (n - 2) * (n - 3)))
    a_kurt <- (6 + (8 / s_kurt) * (2 / s_kurt + sqrt(1 + 4 / (s_kurt^2))))
    z2 <- ((1 - 2 / (9 * a_kurt)) - ((1 - 2 / a_kurt) / (1 + x_kurt * sqrt(2 / (a_kurt - 4))))^(1/3)) / sqrt(2 / (9 * a_kurt))
    
    k2 <- z1^2 + z2^2
    zc <- -qnorm(exp(-0.5 * k2))
    zt <- 0.55 * (n^0.2) - 0.21
    a1 <- (-5 + 3.46 * log(n)) * exp(-1.37 * log(n))
    b1 <- 1 + (0.854 - 0.148 * log(n)) * exp(-0.55 * log(n))
    a2 <- a1 - (2.13 / (1 - 2.37 * log(n))) * zt
    b2 <- 2.13 / (1 - 2.37 * log(n)) + b1
    
    if (zc < -1) z_adj <- zc else if (zc < zt) z_adj <- a1 + b1 * zc else z_adj <- a2 + b2 * zc
    p_val <- 1 - pnorm(z_adj)
    final_chi2 <- -2 * log(p_val)
    return(c(chi2 = final_chi2, p = p_val))
  }
  
  # B. Método Jarque-Bera (Asimetría y Curtosis asintótica estándar)
  jarque_bera_test <- function(v) {
    n <- length(v)
    m1 <- mean(v)
    m2 <- sum((v - m1)^2) / n
    if (m2 == 0) return(c(chi2 = NA, p = NA))
    m3 <- sum((v - m1)^3) / n
    m4 <- sum((v - m1)^4) / n
    
    skew <- m3 / (m2^(3/2))
    kurt <- m4 / (m2^2)
    jb_stat <- (n / 6) * (skew^2 + ((kurt - 3)^2) / 4)
    p_val <- pchisq(jb_stat, df = 2, lower.tail = FALSE)
    return(c(chi2 = jb_stat, p = p_val))
  }
  
  trans_names <- c("cubic", "square", "identity", "square root", "log", 
                   "1/(square root)", "inverse", "1/square", "1/cubic")
  trans_forms <- c("v^3", "v^2", "v", "sqrt(v)", "log(v)", 
                   "1/sqrt(v)", "1/v", "1/(v^2)", "1/(v^3)")
  
  apply_trans <- function(v, type) {
    suppressWarnings(
      switch(type,
             "cubic" = v^3, "square" = v^2, "identity" = v,
             "square root" = sqrt(v), "log" = log(v), "1/(square root)" = 1/sqrt(v),
             "inverse" = 1/v, "1/square" = 1/v^2, "1/cubic" = 1/v^3)
    )
  }
  
  # Matrices de almacenamiento
  res_stata_chi2 <- numeric(9); res_stata_p <- numeric(9)
  res_jb_chi2 <- numeric(9); res_jb_p <- numeric(9)
  res_w_stat <- numeric(9); res_w_p <- numeric(9)
  
  # Muestra para Shapiro si N > 5000 (limitación de R)
  x_shap <- if (n_total > 5000) sample(x, 5000) else x
  
  for (i in 1:9) {
    xt_full <- apply_trans(x, trans_names[i])
    xt_full <- xt_full[is.finite(xt_full) & !is.na(xt_full)]
    
    xt_s <- apply_trans(x_shap, trans_names[i])
    xt_s <- xt_s[is.finite(xt_s) & !is.na(xt_s)]
    
    if (length(xt_full) < 8 || var(xt_full) == 0) {
      res_stata_chi2[i] <- NA; res_stata_p[i] <- NA
      res_jb_chi2[i] <- NA; res_jb_p[i] <- NA
      res_w_stat[i] <- NA; res_w_p[i] <- NA
    } else {
      # 1. Stata sktest
      st_res <- stata_sktest(xt_full)
      res_stata_chi2[i] <- st_res["chi2"]
      res_stata_p[i] <- st_res["p"]
      
      # 2. Jarque-Bera
      jb_res <- jarque_bera_test(xt_full)
      res_jb_chi2[i] <- jb_res["chi2"]
      res_jb_p[i] <- jb_res["p"]
      
      # 3. Shapiro-Wilk
      if (length(xt_s) >= 3) {
        sw_res <- shapiro.test(xt_s)
        res_w_stat[i] <- sw_res$statistic
        res_w_p[i] <- sw_res$p.value
      }
    }
  }
  
  # -------------------------------------------------------------
  # 2. IMPRESIÓN EN CONSOLA CON ENCABEZADOS DE SECCIÓN
  # -------------------------------------------------------------
  cat("\n============================================================\n")
  cat("          LADDER OF POWERS - NORMALITY TESTS COMPARISON     \n")
  cat("============================================================\n")
  
  # SECCIÓN 1: STATA SKTEST
  cat("\n--- 1. AGOSTINO TEST (Chi2 & P-value) ---\n")
  cat(sprintf("%-18s | %-15s | %-12s | %s\n", "Transformation", "formula", "chi2(2)", "P(chi2)"))
  cat(strrep("-", 62), "\n")
  for (i in 1:9) {
    if (is.na(res_stata_chi2[i])) {
      cat(sprintf("%-18s | %-15s | %-12s | %s\n", trans_names[i], trans_forms[i], ".", "."))
    } else {
      c_str <- if (res_stata_chi2[i] > 99999) "." else sprintf("%.2f", res_stata_chi2[i])
      p_str <- if (res_stata_p[i] < 0.0005) "0.000" else sprintf("%.3f", res_stata_p[i])
      cat(sprintf("%-18s | %-15s | %-12s | %s\n", trans_names[i], trans_forms[i], c_str, p_str))
    }
  }
  
  # SECCIÓN 2: SHAPIRO-WILK
  cat("\n--- 2. SHAPIRO-WILK TEST (W & P-value) ---\n")
  cat(sprintf("%-18s | %-15s | %-12s | %s\n", "Transformation", "formula", "W (Stat)", "P(W)"))
  cat(strrep("-", 62), "\n")
  for (i in 1:9) {
    if (is.na(res_w_stat[i])) {
      cat(sprintf("%-18s | %-15s | %-12s | %s\n", trans_names[i], trans_forms[i], ".", "."))
    } else {
      cat(sprintf("%-18s | %-15s | %-12.4f | %.3f\n", trans_names[i], trans_forms[i], res_w_stat[i], res_w_p[i]))
    }
  }
  
  # SECCIÓN 3: JARQUE-BERA
  cat("\n--- 3. JARQUE-BERA TEST (Chi2 & P-value) ---\n")
  cat(sprintf("%-18s | %-15s | %-12s | %s\n", "Transformation", "formula", "JB chi2(2)", "P(JB)"))
  cat(strrep("-", 62), "\n")
  for (i in 1:9) {
    if (is.na(res_jb_chi2[i])) {
      cat(sprintf("%-18s | %-15s | %-12s | %s\n", trans_names[i], trans_forms[i], ".", "."))
    } else {
      cat(sprintf("%-18s | %-15s | %-12.2f | %.3f\n", trans_names[i], trans_forms[i], res_jb_chi2[i], res_jb_p[i]))
    }
  }
  cat("============================================================\n\n")
  
  # Retorno invisible de resultados estilo Stata r()
  r_names <- c("cube", "square", "ident", "sqrt", "log", "invsqrt", "inv", "invsq", "invcube")
  ret_list <- list(N = n_total)
  for (i in 1:9) {
    ret_list[[r_names[i]]] <- res_stata_chi2[i]
    ret_list[[paste0("P_", r_names[i])]] <- res_stata_p[i]
  }
  
  return(invisible(ret_list))
}
