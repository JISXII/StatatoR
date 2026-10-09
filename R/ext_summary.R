#' Extended Summary of Regression Models
#'
#' @description Model-specific summary for lm, binomial glm, nnet::multinom,
#' MASS::polr (logistic/probit), survival::clogit and mlogit::mlogit (ordinary
#' or nested logit), ivreg 2SLS, Poisson glm, MASS::glm.nb (NB2), and fixest
#' linear/IV, Poisson, binomial logit/probit and negative-binomial models.
#' Estimation options belong to the fitting function.
#' @param model Fitted model. Store model frames and Hessians when fitting.
#' @param x An ext_summary result to print.
#' @param ... Arguments forwarded by S3 methods to the summary calculation;
#' additional print arguments are ignored.
#' @param vce "standard", "robust", or "cluster".
#' @param eform Exponentiate logistic slope coefficients; auxiliary parameters
#' remain on their original scale. Multinomial output reports relative risk ratios.
#' @param level Confidence level, proportion or percentage.
#' @param digits Decimal places for display only, from 0 to 12.
#' @param vcov Named covariance matrix on the original parameter scale.
#' @param hc_type Optional HC variant for lm, binomial/Poisson glm or fixest OLS. Default HC1 for
#' lm; HC0 with N/(N-1) correction for unweighted binary glm.
#' @param correction Optional positive multiplier for robust/cluster covariance.
#' Defaults: N/(N-1) for independent discrete outcomes, G/(G-1) for clusters
#' of discrete outcomes, and G/(G-1)*(N-1)/(N-rank) for clustered lm.
#' It replaces the default multiplier except for fixest, where it multiplies
#' the covariance produced by the selected SSC controls.
#' @param cluster One grouping vector aligned to fitted observations (choice
#' cases for mlogit). Clogit defaults to its strata for vce="robust".
#' @param df Inference degrees of freedom: positive finite number for t/F or
#' Inf for z/chi2. Defaults: residual df for lm, G-1 for clustered lm, Inf
#' for discrete models. Custom matrices do not imply cluster df automatically.
#' @param test "auto", "LR", "Wald", or "none". Auto uses F for lm, Wald
#' for alternative-specific/nested logit and robust/custom inference, LR otherwise.
#' @param null_model Optional fitted null model, with the same sample, outcome,
#' likelihood and nested specification. The main model is never refitted.
#' NB2 re-estimates an intercept/offset-only null to maximize its dispersion.
#' @param test_terms Optional exact parameter names for the global Wald test.
#' @param bic_n Optional sample size for BIC. Defaults to rows for clogit/ASCLOGIT
#' (as in the lecture/Stata) and native observation counts for other models.
#' @param vcov_label Optional display label for the covariance estimator.
#' @param show_auxiliary Display ordinal cutpoints/nested dissimilarity parameters.
#' @param auxiliary_tests Test auxiliary parameters against zero; default FALSE.
#' @param vcov_method "auto", "model", "observed", or "opg". Auto uses an
#' observed information matrix for polr and nested logit, native covariance
#' otherwise. observed computes exact ordinal information, joint NB2 information, or numerically
#' differentiates nested logit scores; requires numDeriv for nested models.
#' @param hessian_control Named method.args list for numDeriv::jacobian in
#' nested observed information (for example list(r=6)).
#' @param iia_model Optional ordinary mlogit with identical utility specification
#' and estimation sample, for the nested LR test of all dissimilarities = 1.
#' @param show_deviance Display native null/residual deviance and their degrees
#' of freedom for supported glm models. Default FALSE; requires a glm model.
#' @param only_test Print only the global LR/Wald/F test. The returned object
#' still contains the complete summary. Default FALSE.
#' @param small Use Stata ivregress small-sample covariance and t/F inference
#' for 2SLS. Default FALSE (uncorrected covariance and z/Wald chi-square).
#' @param diagnostics Include conventional ivreg diagnostic tests. Requires
#' ivreg and standard covariance; these are not robust weak-instrument tests.
#' @param fixest_ssc Optional fixest::ssc object controlling parameter counts,
#' cluster corrections and reference degrees of freedom for fixest models.
#' @param fixest_vcov Optional fixest covariance request (string, formula or
#' request object). For example ~id+year or fixest::vcov_NW(). Mutually
#' exclusive with vce/cluster/hc_type/correction/custom vcov options.
#' @param poisson_model Optional fitted Poisson with the same sample, regressors,
#' offsets and fixed effects as an NB2 model. Reports the boundary LR test of
#' alpha=0 using the 50:50 chi-square(0)/chi-square(1) mixture. Conventional only.
#' @return Invisibly, ext_summary with coefficients, vcov, global_test, fit,
#' original model, and display metadata. Numbers are stored without rounding.
#' @details For 2SLS, default inference follows ivregress 2sls: RSS/N model
#' variance, uncorrected robust/cluster scores, z and Wald chi-square. small=TRUE
#' applies N/(N-K) independently and G/(G-1)*(N-1)/(N-K) for one cluster
#' dimension (verified against Stata 19.5), with t/F reference distributions.
#' IV likelihood measures are deliberately unavailable. NB2 uses alpha=1/theta
#' and joint observed information for beta and lnalpha. Count models currently
#' require unweighted integer outcomes. fixest supports ordinary slopes and
#' absorbed categorical effects, not varying slopes or custom nonlinear models.
#' Explicit fixest_ssc delegates correction choices to fixest; these may differ
#' across Stata commands (regress, areg, xtreg, reghdfe, ivregress).
#' @examples
#' m <- lm(mpg ~ wt + am, data=mtcars)
#' ext_summary(m)
#' ext_summary(m, vce="robust", only_test=TRUE)
#' p <- glm(carb ~ wt, data=mtcars, family=poisson())
#' ext_summary(p, eform=TRUE)
#' if (requireNamespace("MASS", quietly=TRUE)) {
#'   nb <- MASS::glm.nb(Days ~ Sex + Age, data=MASS::quine)
#'   ext_summary(nb)
#' }
#' if (requireNamespace("fixest", quietly=TRUE)) {
#'   fe <- fixest::feols(mpg ~ wt | cyl, data=mtcars)
#'   ext_summary(fe, vce="cluster", cluster=~cyl)
#' }
#' @export
ext_summary <- function(model, vce = c("standard", "robust", "cluster"),
                        eform = FALSE, level = .95, digits = 6, vcov = NULL,
                        hc_type = NULL, correction = NULL, cluster = NULL,
                        df = NULL, test = c("auto", "LR", "Wald", "none"),
                        null_model = NULL, test_terms = NULL, bic_n = NULL,
                        vcov_label = NULL, show_auxiliary = TRUE,
                        auxiliary_tests = FALSE,
                        vcov_method = c("auto", "model", "observed", "opg"),
                        hessian_control = list(), iia_model = NULL,
                        show_deviance = FALSE, only_test = FALSE,
                        small = FALSE, diagnostics = FALSE,
                        fixest_ssc = NULL, fixest_vcov = NULL,
                        poisson_model = NULL) {
  UseMethod("ext_summary")
}

#' @rdname ext_summary
#' @export
ext_summary.lm <- function(model, ...) {
  if (inherits(model, "mlm")) stop("Multivariate linear models are not supported.", call. = FALSE)
  .ext_summarize(model, .ext_lm(model), ...)
}
#' @rdname ext_summary
#' @export
ext_summary.glm <- function(model, ...) {
  fam <- stats::family(model)
  if (!((fam$family == "binomial" && fam$link %in% c("logit", "probit")) ||
        (fam$family == "poisson" && fam$link == "log")))
    stop("Supported glm models: binomial logit/probit and Poisson with log link.", call. = FALSE)
  .ext_summarize(model, .ext_glm(model), ...)
}
#' @rdname ext_summary
#' @export
ext_summary.multinom <- function(model, ...) {
  requireNamespace("nnet", quietly = TRUE)
  .ext_summarize(model, .ext_multinom(model), ...)
}
#' @rdname ext_summary
#' @export
ext_summary.polr <- function(model, ...) {
  requireNamespace("MASS", quietly = TRUE)
  if (!model$method %in% c("logistic", "probit"))
    stop("Supported polr links: logistic and probit.", call. = FALSE)
  .ext_summarize(model, .ext_polr(model), ...)
}
#' @rdname ext_summary
#' @export
ext_summary.clogit <- function(model, ...) {
  requireNamespace("survival", quietly = TRUE)
  .ext_summarize(model, .ext_clogit(model), ...)
}
#' @rdname ext_summary
#' @export
ext_summary.mlogit <- function(model, ...) {
  requireNamespace("mlogit", quietly = TRUE)
  .ext_summarize(model, .ext_mlogit(model), ...)
}
#' @rdname ext_summary
#' @export
ext_summary.ivreg <- function(model, ...) {
  .ext_summarize(model, .ext_ivreg(model), ...)
}
#' @rdname ext_summary
#' @export
ext_summary.negbin <- function(model, ...) {
  .ext_summarize(model, .ext_negbin(model), ...)
}
#' @rdname ext_summary
#' @export
ext_summary.fixest <- function(model, ...) {
  .ext_summarize(model, .ext_fixest(model), ...)
}
#' @rdname ext_summary
#' @export
ext_summary.default <- function(model, ...) {
  stop("Supported models: lm, binomial/Poisson glm, multinom, polr, clogit, mlogit, ivreg 2SLS, negbin and supported fixest models.", call. = FALSE)
}

.ext_descriptor <- function(type, b, n, ll, k, ll0 = NA_real_, null_df = 0,
                            auxiliary = character(), global = names(b),
                            eform = TRUE, row_ids = NULL, notes = character()) {
  list(type = type, b = b, n = n, ll = ll, k = k, ll0 = ll0,
       null_df = null_df, auxiliary = auxiliary, global = global,
       eform = eform, row_ids = row_ids, notes = notes,
       group = rep("", length(b)), term = names(b))
}
.ext_lm <- function(m) {
  b <- stats::coef(m); ll <- stats::logLik(m)
  d <- .ext_descriptor("linear", b, stats::nobs(m), as.numeric(ll), attr(ll, "df"),
                       global = setdiff(names(b)[!is.na(b)], "(Intercept)"), eform = FALSE)
  d$residual_df <- stats::df.residual(m)
  d$sample <- stats::model.frame(m)
  d$row_ids <- rownames(d$sample)
  d
}
.ext_glm <- function(m) {
  b <- stats::coef(m); ll <- stats::logLik(m)
  type <- if (m$family$family == "poisson") "poisson" else m$family$link
  d <- .ext_descriptor(type, b, stats::nobs(m), as.numeric(ll), attr(ll, "df"),
                       as.numeric(ll) - (m$null.deviance - m$deviance)/2,
                       stats::nobs(m) - m$df.null,
                       global = setdiff(names(b)[!is.na(b)], "(Intercept)"),
                       eform = m$family$link %in% c("logit", "log"))
  d$sample <- stats::model.frame(m); d$row_ids <- rownames(d$sample)
  d$residual_df <- stats::df.residual(m)
  if (type == "poisson") {
    y <- stats::model.response(d$sample)
    if (any(y < 0) || any(abs(y-round(y)) > 1e-8) || any(m$prior.weights != 1))
      stop("Poisson currently requires unweighted nonnegative integer counts.", call. = FALSE)
  }
  if (!isTRUE(m$converged)) d$notes <- "Model did not converge; inference may be unreliable."
  d
}
# Unweighted 2SLS. Use structural residuals, not second-stage fitted-X residuals.
.ext_ivreg <- function(m) {
  if (!is.null(m$method) && m$method != "OLS")
    stop("ivreg support is 2SLS (method='OLS'), not 2SM/2SMM/LIML/GMM.", call. = FALSE)
  mf <- stats::model.frame(m); w <- stats::model.weights(mf)
  if (!is.null(w) && any(w != 1)) stop("Weighted IV is not yet supported.", call. = FALSE)
  if (!is.null(stats::model.offset(mf))) stop("IV offsets are not yet supported.", call. = FALSE)
  b <- stats::coef(m); nm <- names(b)[!is.na(b)]
  X <- stats::model.matrix(m, component = "regressors")[, nm, drop = FALSE]
  Z <- stats::model.matrix(m, component = "instruments")
  qz <- qr(Z); Xhat <- qr.fitted(qz, X)
  if (qr(Xhat)$rank < ncol(Xhat)) stop("Instruments do not identify all retained regressors.", call. = FALSE)
  y <- as.numeric(stats::model.response(mf)); u <- as.numeric(y-X %*% b[nm])
  A <- solve(crossprod(Xhat)); dimnames(A) <- list(nm, nm)
  cte <- as.integer("(Intercept)" %in% nm); n <- nrow(X); k <- ncol(X)
  tss <- sum((y-if (cte) mean(y) else 0)^2); rss <- sum(u^2)
  d <- .ext_descriptor("iv_2sls", b, n, NA_real_, k, eform = FALSE,
                       global = setdiff(nm, "(Intercept)"), row_ids = rownames(mf),
                       notes = "2SLS uses structural residuals. LR, likelihood AIC/BIC and pseudo-R2 are not reported.")
  d$sample <- mf; d$residual_df <- n-k; d$force_wald <- TRUE
  d$bread <- A; d$scores <- Xhat*u; d$rss <- rss
  d$fit_extra <- list(r_squared = if (tss > 0) 1-rss/tss else NA_real_,
                     adj_r_squared = if (tss > 0) 1-rss/tss*(n-cte)/(n-k) else NA_real_,
                     residual_df = n-k)
  d$extension <- TRUE; d
}

.ext_negbin <- function(m) {
  if (m$family$link != "log") stop("NB2 support requires a log link.", call. = FALSE)
  mf <- stats::model.frame(m); y <- stats::model.response(mf); w <- stats::model.weights(mf)
  if ((!is.null(w) && any(w != 1)) || any(y < 0) || any(abs(y-round(y)) > 1e-8))
    stop("NB2 currently requires unweighted nonnegative integer counts.", call. = FALSE)
  if (!is.finite(m$theta) || m$theta <= 0) stop("A positive finite NB2 theta is required.", call. = FALSE)
  if ("lnalpha" %in% names(stats::coef(m))) stop("lnalpha is reserved for the NB2 dispersion parameter.", call. = FALSE)
  b <- c(stats::coef(m), lnalpha = -log(m$theta)); ll <- stats::logLik(m)
  d <- .ext_descriptor("negbin", b, nrow(mf), as.numeric(ll), sum(!is.na(b)),
                       auxiliary = "lnalpha", global = setdiff(names(b)[!is.na(b)], c("(Intercept)", "lnalpha")),
                       row_ids = rownames(mf))
  d$sample <- mf; d$residual_df <- m$df.residual
  # glm.nb null.deviance holds theta fixed: it is NOT the maximized NB2 null LL.
  # Re-estimate dispersion in an intercept/offset-only model on the exact sample.
  has_intercept <- attr(m$terms, "intercept") == 1
  off <- stats::model.offset(mf); if (is.null(off)) off <- rep(0, length(y))
  null_data <- data.frame(y = y, off = off)
  null <- tryCatch(suppressWarnings(MASS::glm.nb(
    if (has_intercept) y~1+offset(off) else y~0+offset(off),
    data = null_data, control = stats::glm.control(epsilon=1e-10, maxit=200))), error=function(e) NULL)
  if (!is.null(null) && isTRUE(null$converged) && is.null(null$th.warn)) {
    d$ll0 <- as.numeric(stats::logLik(null)); d$null_df <- has_intercept+1L
  } else {
    d$force_wald <- TRUE
    d$notes <- "Maximized NB2 null likelihood unavailable; supply null_model for an LR comparison."
  }
  if (!isTRUE(m$converged) || !is.null(m$th.warn))
    d$notes <- c(d$notes, "NB2 estimation/dispersion did not converge; inspect the fitted model.")
  d
}

.ext_negbin_information <- function(m) {
  if (!requireNamespace("numDeriv", quietly=TRUE)) stop("Joint NB2 observed information requires numDeriv.", call.=FALSE)
  beta <- stats::coef(m); beta <- beta[!is.na(beta)]
  X <- stats::model.matrix(m)[, names(beta), drop=FALSE]
  y <- as.numeric(stats::model.response(stats::model.frame(m)))
  off <- stats::model.offset(stats::model.frame(m)); if (is.null(off)) off <- rep(0,length(y))
  b <- c(beta, lnalpha=-log(m$theta)); k <- length(beta)
  scores <- function(par) {
    mu <- exp(as.numeric(X %*% par[seq_len(k)])+off); theta <- exp(-par[k+1L])
    sb <- X*as.numeric((y-mu)/(1+mu/theta))
    sa <- theta*(digamma(theta)-digamma(y+theta)+log1p(mu/theta)+(y-mu)/(theta+mu))
    S <- cbind(sb, lnalpha=sa); colnames(S) <- names(b); S
  }
  S <- scores(b)
  H <- -numDeriv::jacobian(function(par) colSums(scores(par)), b)
  H <- (H+t(H))/2; dimnames(H) <- list(names(b),names(b))
  if (any(!is.finite(H)) || min(eigen(H,symmetric=TRUE,only.values=TRUE)$values) <= 0)
    stop("NB2 observed information is not positive definite; inspect dispersion/identification or supply vcov.",call.=FALSE)
  list(vcov=solve(H),scores=S,information=H)
}

.ext_fixest <- function(m) {
  if (!requireNamespace("fixest", quietly=TRUE)) stop("fixest is required.",call.=FALSE)
  if (!is.null(m$fixef_terms) || !is.null(m$slope_flag))
    stop("Varying-slope fixed effects are not yet supported.",call.=FALSE)
  w <- m$weights
  if (!is.null(w) && any(w != 1)) stop("Weighted fixest models are not yet supported.",call.=FALSE)
  iv <- isTRUE(m$is_iv)
  if (iv && !identical(m$iv_stage,2)) stop("Use the second-stage fixest IV model.",call.=FALSE)
  b <- stats::coef(m); family <- m$family
  if (m$method_type == "feols") type <- if (iv) "fixest_iv" else "fixest_linear"
  else if (inherits(family,"family") && family$family == "poisson" && family$link == "log") type <- "fixest_poisson"
  else if (inherits(family,"family") && family$family == "binomial" && family$link %in% c("logit","probit"))
    type <- paste0("fixest_",family$link)
  else if (is.character(family) && family == "negbin") type <- "fixest_negbin"
  else if (is.character(family) && family %in% c("poisson","logit")) type <- paste0("fixest_",family)
  else stop("fixest support: feols OLS/2SLS, Poisson, binomial logit/probit, and NB2; nonlinear/custom likelihoods are not supported.",call.=FALSE)
  n <- stats::nobs(m); ll <- if (iv) NA_real_ else as.numeric(stats::logLik(m))
  ids <- as.character(fixest::obs(m))
  k <- m$nparams + as.integer(type == "fixest_negbin") + as.integer(type == "fixest_linear")
  if (type == "fixest_negbin") b <- c(b,lnalpha=-log(as.numeric(m$theta)))
  d <- .ext_descriptor(type,b,n,ll,k,
                       ll0=if (iv || is.null(m$ll_null)) NA_real_ else m$ll_null,
                       auxiliary=if (type == "fixest_negbin") "lnalpha" else character(),
                       global=setdiff(names(b),c("(Intercept)","lnalpha")),
                       eform=!type %in% c("fixest_linear","fixest_iv","fixest_probit"),row_ids=ids)
  d$extension <- TRUE; d$force_wald <- TRUE; d$residual_df <- n-m$nparams
  if (type %in% c("fixest_poisson","fixest_negbin")) {
    y <- as.numeric(stats::model.matrix(m,type="lhs"))
    if (any(y < 0) || any(abs(y-round(y)) > 1e-8)) stop("Count models require nonnegative integer outcomes.",call.=FALSE)
  }
  if (type %in% c("fixest_logit","fixest_probit")) {
    y <- if (!is.null(m$y)) as.numeric(m$y) else as.numeric(stats::model.matrix(m,type="lhs"))
    if (!all(y %in% c(0,1))) stop("fixest binomial support requires individual binary outcomes.",call.=FALSE)
  }
  d$fit_extra <- list(fixed_effects=m$fixef_sizes)
  if (type %in% c("fixest_linear","fixest_iv")) {
    fs <- fixest::fitstat(m,c("r2","ar2","wr2","war2"))
    d$fit_extra <- c(d$fit_extra,list(r_squared=fs$r2,adj_r_squared=fs$ar2,
                                    within_r_squared=fs$wr2,adj_within_r_squared=fs$war2,
                                    residual_df=d$residual_df))
    d$rss <- sum(m$residuals^2)
  }
  if (length(m$fixef_sizes)) d$notes <- c(d$notes,
    "Absorbed fixed effects are not printed as coefficients; the global test covers the displayed slopes only.")
  if (iv) d$notes <- c(d$notes,"2SLS: LR, likelihood AIC/BIC and pseudo-R2 are not reported.")
  if (!is.null(m$convStatus) && !isTRUE(m$convStatus)) d$notes <- c(d$notes,"Model did not converge.")
  d
}

.ext_extension_covariance <- function(m,d,vce,vcov,hc_type,correction,cluster,small,df,
                                      fixest_ssc,fixest_vcov,vcov_method) {
  iv <- d$type %in% c("iv_2sls","fixest_iv"); linear <- d$type == "fixest_linear"
  G <- NULL; meta <- NULL
  if (vcov_method != "auto" && vcov_method != "model")
    stop("IV/fixest use model covariance or their covariance controls; observed/opg are not supported here.",call.=FALSE)
  if (!is.null(vcov)) {
    V <- vcov; label <- "Custom"; infer_df <- if (iv) if (small) d$residual_df else Inf else if (linear) d$residual_df else Inf
  } else if (d$type == "iv_2sls") {
    if (!is.null(hc_type)) stop("2SLS uses projected-regressor scores; hc_type is not supported. Use small/correction or explicit vcov.",call.=FALSE)
    n <- d$n; k <- sum(!is.na(d$b)); factor <- if (small) n/(n-k) else 1
    if (vce == "standard") { V <- d$bread*d$rss/n*factor; label <- if (small) "Standard (small)" else "Standard (asymptotic)" }
    else {
      S <- d$scores
      if (!is.null(cluster)) {
        cluster <- .ext_cluster_vector(cluster,d$row_ids); G <- length(unique(cluster))
        S <- rowsum(S,cluster,reorder=FALSE)
        factor <- if (small) (n-1)/(n-k)*G/(G-1) else 1
      } else if (vce == "cluster") stop("Supply cluster aligned to the IV estimation sample.",call.=FALSE)
      if (!is.null(correction)) factor <- correction
      V <- d$bread %*% crossprod(S) %*% d$bread*factor
      label <- if (is.null(G)) "Robust (2SLS scores)" else paste0("Cluster (2SLS, G=",G,")")
    }
    infer_df <- if (small) if (is.null(G)) d$residual_df else G-1 else Inf
  } else {
    if (!is.null(fixest_vcov) && (vce != "standard" || !is.null(cluster) || !is.null(hc_type) || !is.null(correction)))
      stop("Supply fixest_vcov without vce/cluster/hc_type/correction.",call.=FALSE)
    if (!is.null(hc_type) && (!linear || vce != "robust" || !is.null(cluster)))
      stop("fixest hc_type is available for independent robust linear OLS only.",call.=FALSE)
    request <- if (!is.null(fixest_vcov)) fixest_vcov else if (!is.null(cluster)) "cluster" else
      switch(vce,standard="iid",robust=if (is.null(hc_type)) "hetero" else match.arg(hc_type,c("HC1","HC2","HC3")),cluster="cluster")
    if (is.null(fixest_ssc)) {
      # Explicit settings avoid global setFixest_ssc() silently changing this summary.
      ssc <- fixest::ssc(K.adj=linear, K.fixef="nonnested",K.exact=TRUE,
                         G.adj=!iv,G.df="min",t.df="min")
    } else {
      if (!inherits(fixest_ssc,"ssc_type")) stop("fixest_ssc must be created with fixest::ssc().",call.=FALSE)
      ssc <- fixest_ssc
    }
    work_model <- m
    # Native NB2 scores need the nuisance-FE projection for robust inference.
    # Replace scores on a local copy; coefficients/fit remain untouched.
    if (d$type == "fixest_negbin" && length(m$fixef_id))
      work_model$scores <- .ext_fixest_nb_scores(m)
    args <- list(object=work_model,vcov=request,ssc=ssc,attr=TRUE,vcov_fix=FALSE)
    if (!is.null(cluster)) {
      if (inherits(cluster,"formula")) args$vcov <- cluster else {
        cluster <- .ext_cluster_vector(cluster,d$row_ids)
        args$vcov <- NULL; args$cluster <- list(cluster)
      }
    }
    raw <- do.call(stats::vcov,args); meta <- attributes(raw)
    if (!is.null(meta$G)) G <- meta$G
    else if (!is.null(meta$min_cluster_size)) G <- meta$min_cluster_size
    else if (!is.null(cluster) && !inherits(cluster,"formula")) G <- length(unique(cluster))
    V <- raw
    if (d$type == "fixest_negbin") {
      nm <- names(d$b); orig <- c(setdiff(nm,"lnalpha"),".theta")
      J <- diag(c(rep(1,length(nm)-1L),-1/as.numeric(m$theta)))
      V <- J %*% raw[orig,orig,drop=FALSE] %*% J; dimnames(V) <- list(nm,nm)
    }
    # fixest IID linear covariance always uses residual df; ivregress defaults to RSS/N.
    if (iv && is.null(fixest_ssc)) {
      if (small && !is.null(G) && grepl(" & ",meta$vcov_type,fixed=TRUE))
        stop("IV small-sample multiway clustering requires explicit fixest_ssc; automatic ivregress small correction covers one cluster dimension.",call.=FALSE)
      iid <- identical(meta$vcov_type,"IID")
      if (iid) {
        A <- solve(m$hessian)
        dimnames(A) <- list(names(stats::coef(m)),names(stats::coef(m)))
        V <- A*d$rss/if (small) d$residual_df else d$n
      }
      if (!iid && small) V <- V*if (is.null(G)) d$n/d$residual_df else (d$n-1)/d$residual_df*G/(G-1)
    }
    if (!linear && !iv && is.null(fixest_ssc) && is.null(fixest_vcov) &&
        vce == "robust" && is.null(cluster)) V <- V*d$n/(d$n-1)
    if (!is.null(correction)) V <- V*correction
    label <- paste0("fixest: ",meta$vcov_type)
    if (!is.null(fixest_ssc)) label <- paste0(label," (custom SSC)")
    infer_df <- if (iv) if (small) if (is.null(G)) d$residual_df else G-1 else Inf else
      if (linear) meta$df.t else Inf
  }
  if (!is.null(df)) infer_df <- df
  list(vcov=V,label=label,df=infer_df,clusters=G,metadata=meta,
       altered=!is.null(vcov) || vce != "standard" || !is.null(cluster) || !is.null(fixest_vcov))
}

.ext_cluster_vector <- function(cluster,ids) {
  if (!is.atomic(cluster) || anyNA(cluster) || length(cluster) != length(ids))
    stop("cluster must be one vector aligned to the estimation sample, without missing values.",call.=FALSE)
  if (!is.null(names(cluster))) {
    if (anyDuplicated(names(cluster)) || !setequal(names(cluster),ids)) stop("Named cluster vector must match estimation row identifiers.",call.=FALSE)
    cluster <- cluster[ids]
  }
  if (length(unique(cluster)) < 2) stop("At least two clusters are required.",call.=FALSE)
  cluster
}

.ext_fixest_nb_scores <- function(m) {
  X <- stats::model.matrix(m)[,names(stats::coef(m)),drop=FALSE]
  y <- as.numeric(stats::model.matrix(m,type="lhs")); mu <- as.numeric(stats::fitted(m))
  theta <- as.numeric(m$theta); alpha <- 1/theta
  q <- (y-mu)/(1+alpha*mu)
  w <- mu*theta*(theta+y)/(theta+mu)^2
  cross <- (y-mu)*alpha*mu/(1+alpha*mu)^2
  z <- cross/w
  centered <- fixest::demean(cbind(X,z),f=m$fixef_id,weights=w,tol=1e-10)
  s_alpha <- theta*(digamma(theta)-digamma(y+theta)+log1p(mu/theta)+(y-mu)/(theta+mu))
  # Efficient per-observation scores: S_a - S_FE H_FE,FE^-1 H_FE,a.
  S <- cbind(centered[,seq_len(ncol(X)),drop=FALSE]*q,
             (s_alpha-q*(z-centered[,ncol(X)+1L]))*(-1/theta))
  colnames(S) <- c(colnames(X),".theta"); S
}

.ext_multinom <- function(m) {
  if (isTRUE(m$censored) || (!is.null(m$call$summ) && as.character(m$call$summ) != "0") ||
      (!is.null(m$decay) && m$decay != 0))
    stop("multinom support requires uncensored outcomes, summ=0 and decay=0.", call. = FALSE)
  if (is.null(m$model)) stop("Fit multinom with model=TRUE to preserve its estimation sample.", call. = FALSE)
  mf <- m$model; y <- stats::model.response(mf)
  if (!is.factor(y)) stop("multinom currently requires an individual factor response (not counts).", call. = FALSE)
  if (is.null(m$Hessian)) stop("Fit multinom with Hess=TRUE.", call. = FALSE)
  B <- stats::coef(m)
  binary <- is.null(dim(B))
  if (is.null(dim(B))) B <- matrix(B, 1, dimnames = list(m$lev[2], names(B)))
  nm <- as.vector(t(outer(rownames(B), colnames(B), paste, sep = ":")))
  if (binary) nm <- colnames(B)
  b <- stats::setNames(as.vector(t(B)), nm)
  weights <- as.numeric(m$weights)
  if (!length(weights)) weights <- rep(1, nrow(mf))
  counts <- tapply(weights, y, sum)
  intercept <- attr(m$terms, "intercept")
  ll0 <- if (!is.null(stats::model.offset(mf))) NA_real_ else
    if (intercept) sum(counts[counts > 0] * log(counts[counts > 0]/sum(counts)))
    else -sum(counts) * log(length(m$lev))
  d <- .ext_descriptor("multinomial", b, nrow(mf), as.numeric(stats::logLik(m)), m$edf,
                       ll0, if (intercept) length(m$lev)-1L else 0,
                       global = nm[rep(colnames(B) != "(Intercept)", nrow(B))],
                       row_ids = rownames(mf))
  d$group <- rep(rownames(B), each = ncol(B)); d$term <- rep(colnames(B), nrow(B))
  d$base <- m$lev[1]; d$sample <- mf; d$weights <- weights
  if (m$convergence != 0) d$notes <- "Model did not converge; inference may be unreliable."
  d
}
.ext_polr <- function(m) {
  if (is.null(m$model) || is.null(m$Hessian)) stop("Fit polr with model=TRUE and Hess=TRUE.", call. = FALSE)
  b <- c(stats::coef(m), m$zeta); mf <- m$model
  w <- stats::model.weights(mf)
  if (is.null(w)) w <- rep(1, nrow(mf))
  y <- stats::model.response(mf); counts <- tapply(w, y, sum)
  ll0 <- if (is.null(stats::model.offset(mf)))
    sum(counts[counts > 0] * log(counts[counts > 0]/sum(counts))) else NA_real_
  d <- .ext_descriptor(if (m$method == "logistic") "ordinal_logit" else "ordinal_probit",
                       b, stats::nobs(m), as.numeric(stats::logLik(m)), m$edf,
                       ll0, length(m$zeta), auxiliary = names(m$zeta),
                       global = names(stats::coef(m)), eform = m$method == "logistic",
                       row_ids = rownames(mf))
  d$sample <- mf; d$weights <- w
  d$group <- c(rep("Slopes", length(m$coefficients)), rep("Cutpoints", length(m$zeta)))
  d$term <- c(names(m$coefficients), paste0("/cut", seq_along(m$zeta)))
  if (m$convergence != 0) d$notes <- "Model did not converge; inference may be unreliable."
  d
}
.ext_clogit <- function(m) {
  if (is.null(m$model) || is.null(m$x) || is.null(m$y))
    stop("Fit clogit with model=TRUE, x=TRUE and y=TRUE.", call. = FALSE)
  b <- stats::coef(m)
  d <- .ext_descriptor("conditional", b, m$n, m$loglik[2], sum(!is.na(b)),
                       m$loglik[1], row_ids = rownames(m$x))
  d$sample <- m$model; d$groups <- m$strata
  if (is.null(d$groups)) d$groups <- factor(rep(1, m$n))
  sizes <- table(d$groups); events <- rowsum(m$y[, ncol(m$y)], d$groups)
  informative <- events[, 1] > 0 & events[, 1] < as.numeric(sizes)
  if (any(!informative)) stop("Remove noninformative clogit strata before fitting (all outcomes 0 or all 1), so sample counts and corrections match Stata.", call. = FALSE)
  d$fit_extra <- list(cases = length(sizes), informative_cases = sum(informative),
                      alternatives = c(min(sizes), mean(sizes), max(sizes)))
  if (any(!informative)) d$notes <- "Some strata contain no within-stratum outcome variation."
  d
}
.ext_choice_index <- function(m) {
  if (inherits(m$model, "dfidx")) {
    requireNamespace("dfidx", quietly = TRUE)
    ix <- dfidx::idx(m$model)
  } else ix <- attr(m$model, "index")
  if (is.null(ix) || ncol(ix) < 2) stop("Cannot recover choice-case and alternative indexes.", call. = FALSE)
  # dfidx may have a second case index (panel id); the alternative is last.
  data.frame(case = as.character(ix[[1]]), alternative = as.character(ix[[ncol(ix)]]),
             stringsAsFactors = FALSE)
}
.ext_mlogit <- function(m) {
  if (!is.null(m$rpar) || isTRUE(m$call$probit) || isTRUE(m$call$heterosc) ||
      isTRUE(m$call$ranked) || isTRUE(m$call$panel))
    stop("Only ordinary and nested mlogit models are supported; not mixed/probit/heterosc/ranked/panel fits.", call. = FALSE)
  b <- stats::coef(m)
  fixed <- attr(m$coefficients, "fixed")
  fixed_names <- names(fixed)[fixed %in% TRUE]
  ix <- .ext_choice_index(m); cases <- unique(ix$case); sizes <- table(ix$case)
  aux <- grep("^iv([.:]|$)", names(b), value = TRUE)
  is_nested <- !is.null(m$nests)
  if (is_nested && length(fixed_names)) stop("Nested models with fixed parameters require a separate adapter.", call. = FALSE)
  ll <- stats::logLik(m)
  null_ll <- if (inherits(m, "micsr")) unname(m$logLik["null"]) else {
    freq <- as.numeric(m$freq); sum(freq[freq > 0] * log(freq[freq > 0]/sum(freq)))
  }
  if (length(unique(sizes)) != 1) null_ll <- NA_real_
  has_intercepts <- any(grepl("^\\(Intercept\\)(:|$)", names(b)))
  if (!has_intercepts) null_ll <- -sum(log(as.numeric(sizes)))
  pure_multinomial <- !is_nested && all(vapply(names(b), function(nm)
    any(endsWith(nm, paste0(":", unique(ix$alternative)))), logical(1)))
  d <- .ext_descriptor(if (is_nested) "nested" else if (pure_multinomial) "multinomial_choice" else "alternative_specific", b,
                       if (pure_multinomial) length(cases) else nrow(ix), as.numeric(ll), length(b), null_ll,
                       if (has_intercepts) length(unique(ix$alternative))-1L else 0L, auxiliary = aux,
                       global = setdiff(names(b)[!grepl("^\\(Intercept\\)(:|$)", names(b))], aux),
                       eform = !is_nested, row_ids = cases)
  d$sample <- ix; d$ix <- ix
  d$force_wald <- length(fixed_names) > 0
  if (length(fixed_names)) d$notes <- paste("Fixed parameters excluded from inference and parameter count:",
      paste(paste0(fixed_names, "=", m$coefficients[fixed_names]), collapse = ", "))
  d$fit_extra <- list(cases = length(cases), alternatives = c(min(sizes), mean(sizes), max(sizes)))
  nonbase <- sub("^\\(Intercept\\):", "", grep("^\\(Intercept\\):", names(b), value = TRUE))
  represented <- unique(ix$alternative)[vapply(unique(ix$alternative), function(a)
    any(endsWith(names(b)[!names(b) %in% aux], paste0(":", a))), logical(1))]
  d$base <- setdiff(unique(ix$alternative), if (length(nonbase)) nonbase else represented)
  if (length(d$base) != 1L) d$base <- character()
  d$group <- rep("Generic slopes", length(b))
  for (alt in unique(ix$alternative)) {
    rows <- endsWith(names(b), paste0(":", alt)) & !names(b) %in% aux
    d$group[rows] <- alt
    d$term[rows] <- substr(names(b)[rows], 1, nchar(names(b)[rows])-nchar(alt)-1L)
  }
  d$group[names(b) %in% aux] <- "Dissimilarity"
  if (is_nested && any(b[aux] <= 0 | b[aux] > 1))
    d$notes <- "Dissimilarity parameters outside (0, 1]: this fit violates the usual RUM-consistency restriction."
  if (inherits(m, "micsr") && !is.null(m$est.stat$code) && !m$est.stat$code %in% c(1, 2))
    d$notes <- c(d$notes, "Optimizer did not report convergence.")
  d
}

.ext_polr_information <- function(m) {
  X <- stats::model.matrix(m)[, names(m$coefficients), drop = FALSE]
  y <- as.integer(stats::model.response(m$model)); k <- ncol(X); a <- length(m$zeta)
  b <- c(m$coefficients, m$zeta); nm <- names(b)
  lp <- as.numeric(m$lp); cuts <- c(-Inf, m$zeta, Inf)
  w <- stats::model.weights(m$model)
  if (is.null(w)) w <- rep(1, nrow(X))
  info <- matrix(0, k+a, k+a, dimnames = list(nm, nm))
  scores <- matrix(0, nrow(X), k+a, dimnames = list(rownames(X), nm))
  cdf <- if (m$method == "logistic") stats::plogis else stats::pnorm
  density <- if (m$method == "logistic") stats::dlogis else stats::dnorm
  derivative <- if (m$method == "logistic") function(z) density(z)*(1-2*cdf(z)) else
    function(z) ifelse(is.finite(z), -z * density(z), 0)
  for (i in seq_len(nrow(X))) {
    l <- cuts[y[i]] - lp[i]; u <- cuts[y[i]+1L] - lp[i]
    # Compute upper-tail differences in the positive tail to avoid cancellation.
    prob <- if (l > 0) cdf(l, lower.tail = FALSE)-cdf(u, lower.tail = FALSE) else cdf(u)-cdf(l)
    if (prob <= 0 && w[i] > 0) stop("Ordinal probability underflow; inference unavailable.", call. = FALSE)
    if (w[i] == 0) next
    fl <- density(l); fu <- density(u); dl <- derivative(l); du <- derivative(u)
    v <- numeric(k+a); H <- matrix(0, k+a, k+a)
    if (k > 0) {
      v[seq_len(k)] <- (fl-fu)*X[i, ]
      H[seq_len(k), seq_len(k)] <- (du-dl)*tcrossprod(X[i, ])
    }
    if (y[i] > 1) {
      lo <- k+y[i]-1L; v[lo] <- -fl; H[lo, lo] <- -dl
      if (k > 0) H[seq_len(k), lo] <- H[lo, seq_len(k)] <- dl*X[i, ]
    }
    if (y[i] <= a) {
      hi <- k+y[i]; v[hi] <- fu; H[hi, hi] <- du
      if (k > 0) H[seq_len(k), hi] <- H[hi, seq_len(k)] <- -du*X[i, ]
    }
    scores[i, ] <- w[i]*v/prob
    info <- info+w[i]*(tcrossprod(v)/prob^2-H/prob)
  }
  list(vcov = solve(info), scores = scores, information = info)
}

.ext_nested_information <- function(m, control = list()) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Nested observed information requires numDeriv.", call. = FALSE)
  ix <- .ext_choice_index(m); X <- stats::model.matrix(m)
  beta_names <- colnames(X); b <- stats::coef(m); aux <- setdiff(names(b), beta_names)
  if (any(!grepl("^iv([.:]|$)", aux))) stop("Unrecognized nested parameterization.", call. = FALSE)
  nests <- m$nests
  if (is.null(nests)) stop("No nests in this model.", call. = FALSE)
  membership <- vapply(ix$alternative, function(a) {
    hit <- which(vapply(nests, function(z) a %in% z, logical(1)))
    if (length(hit) != 1) stop("Every alternative must belong to exactly one nest.", call. = FALSE)
    hit
  }, integer(1))
  y <- as.logical(stats::model.response(m$model))
  w <- stats::model.weights(m$model)
  if (!is.null(w) && any(w != 1)) stop("Nested observed information currently requires unweighted cases.", call. = FALSE)
  unscaled <- isTRUE(m$call$unscaled)
  case_rows <- split(seq_len(nrow(X)), factor(ix$case, levels = unique(ix$case)))
  if (any(vapply(case_rows, function(i) sum(y[i]) != 1, logical(1)))) stop("One choice per case is required.", call. = FALSE)
  lse <- function(z) { h <- max(z); h+log(sum(exp(z-h))) }
  score <- function(par, details = FALSE) {
    beta <- par[beta_names]; lambda <- par[aux]
    if (length(aux) == 1) lambda <- rep(lambda, length(nests))
    if (length(lambda) != length(nests)) stop("Cannot map nest parameters.", call. = FALSE)
    eta <- as.numeric(X %*% beta)
    S <- matrix(0, length(case_rows), length(par), dimnames = list(names(case_rows), names(par)))
    values <- numeric(length(case_rows))
    for (i in seq_along(case_rows)) {
      rows <- case_rows[[i]]; available <- unique(membership[rows])
      inclusive <- numeric(length(available)); within <- numeric(length(rows))
      xbar <- matrix(0, length(available), length(beta)); etabar <- numeric(length(available))
      for (g in seq_along(available)) {
        a <- available[g]; loc <- which(membership[rows] == a)
        u <- eta[rows[loc]] / if (unscaled) 1 else lambda[a]
        inclusive[g] <- lse(u)
        within[loc] <- exp(u-inclusive[g])
        xbar[g, ] <- colSums(X[rows[loc], , drop = FALSE]*within[loc])
        etabar[g] <- sum(eta[rows[loc]]*within[loc])
      }
      upper <- lambda[available]*inclusive; norm <- lse(upper)
      nest_prob <- exp(upper-norm); chosen <- which(y[rows]); a <- membership[rows[chosen]]
      g <- match(a, available)
      values[i] <- log(within[chosen])+upper[g]-norm
      if (!unscaled) {
        S[i, beta_names] <- X[rows[chosen], ]/lambda[a]+(1-1/lambda[a])*xbar[g, ]-
          colSums(xbar*nest_prob)
        base <- inclusive-etabar/lambda[available]
        d <- -nest_prob*base
        d[g] <- d[g]+inclusive[g]+(etabar[g]-eta[rows[chosen]])/lambda[a]^2-etabar[g]/lambda[a]
      } else {
        S[i, beta_names] <- X[rows[chosen], ]+(lambda[a]-1)*xbar[g, ]-
          colSums(xbar*(nest_prob*lambda[available]))
        d <- -nest_prob*inclusive; d[g] <- d[g]+inclusive[g]
      }
      if (length(aux) == 1) S[i, aux] <- sum(d) else S[i, aux[available]] <- d
    }
    if (details) list(scores = S, ll = sum(values)) else colSums(S)
  }
  initial <- score(b, TRUE)
  if (abs(initial$ll-as.numeric(stats::logLik(m))) > 1e-6)
    stop("Nested likelihood parameterization does not match the fitted model.", call. = FALSE)
  control <- utils::modifyList(list(d = .01, r = 6), control)
  info <- -numDeriv::jacobian(score, b, method = "Richardson", method.args = control)
  info <- (info+t(info))/2; dimnames(info) <- list(names(b), names(b))
  if (min(eigen(info, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("Nested observed information is not positive definite; inspect estimation or supply vcov.", call. = FALSE)
  list(vcov = solve(info), scores = initial$scores, information = info)
}

.ext_scores <- function(m, d) {
  if (d$type == "linear" || d$type %in% c("logit", "probit", "poisson")) {
    if (!requireNamespace("sandwich", quietly = TRUE)) stop("This covariance estimator requires sandwich.", call. = FALSE)
    return(sandwich::estfun(m))
  }
  if (d$type == "multinomial") {
    X <- stats::model.matrix(m$terms, m$model, contrasts.arg = m$contrasts)
    X <- X[, m$vcoefnames, drop = FALSE]
    y <- stats::model.response(m$model); K <- length(m$lev)
    P <- m$fitted.values
    if (is.null(dim(P)) || ncol(P) == 1L) P <- cbind(1-as.numeric(P), as.numeric(P))
    S <- do.call(cbind, lapply(seq.int(2L, K), function(j)
      X * as.numeric(y == m$lev[j])-X*P[, j]))
    S <- S*d$weights; colnames(S) <- names(d$b); rownames(S) <- d$row_ids
    return(S)
  }
  if (d$type %in% c("ordinal_logit", "ordinal_probit")) return(.ext_polr_information(m)$scores)
  if (d$type == "negbin") return(.ext_negbin_information(m)$scores)
  if (d$type == "conditional") {
    if (m$method == "exact") {
      sizes <- table(d$groups); events <- rowsum(m$y[, ncol(m$y)], d$groups)[, 1]
      if (any(events > 1 & events < sizes)) stop("Robust exact clogit with multiple chosen outcomes is not supported; supply vcov.", call. = FALSE)
      # One event per informative stratum: Breslow is the same conditional likelihood.
      m$method <- "breslow"
    }
    return(stats::residuals(m, type = "score"))
  }
  if (d$type == "nested") return(.ext_nested_information(m)$scores)
  if (d$type %in% c("alternative_specific", "multinomial_choice")) {
    S <- m$gradient
    if (is.null(dim(S)) || ncol(S) < length(d$b))
      stop("Per-case scores are required; use mlogit >= 2.0-0 or supply vcov.", call. = FALSE)
    if (is.null(colnames(S))) colnames(S) <- names(m$coefficients)
    return(S[, names(d$b), drop = FALSE])
  }
}
.ext_validate_vcov <- function(V, nm) {
  if (!is.matrix(V) || !is.numeric(V) || is.null(rownames(V)) || is.null(colnames(V)) ||
      anyDuplicated(rownames(V)) || anyDuplicated(colnames(V)) ||
      !all(nm %in% rownames(V)) || !all(nm %in% colnames(V)))
    stop("vcov must be a numeric matrix named for all estimable parameters, including auxiliary parameters.", call. = FALSE)
  V <- V[nm, nm, drop = FALSE]
  V <- matrix(as.numeric(V), nrow(V), dimnames = dimnames(V))
  if (any(!is.finite(V)) || !isTRUE(all.equal(unname(V), unname(t(V)), tolerance = 1e-8)))
    stop("Covariance matrix must be finite and symmetric.", call. = FALSE)
  V <- (V+t(V))/2
  if (any(diag(V) < 0) || min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) <
      -sqrt(.Machine$double.eps)*max(1, max(abs(V))))
    stop("Covariance matrix must be positive semidefinite.", call. = FALSE)
  V
}
.ext_summarize <- function(model, d, vce = c("standard", "robust", "cluster"),
                           eform = FALSE, level = .95, digits = 6, vcov = NULL,
                           hc_type = NULL, correction = NULL, cluster = NULL,
                           df = NULL, test = c("auto", "LR", "Wald", "none"),
                           null_model = NULL, test_terms = NULL, bic_n = NULL,
                           vcov_label = NULL, show_auxiliary = TRUE,
                           auxiliary_tests = FALSE,
                           vcov_method = c("auto", "model", "observed", "opg"),
                           hessian_control = list(), iia_model = NULL,
                           show_deviance = FALSE, only_test = FALSE,
                           small = FALSE, diagnostics = FALSE,
                           fixest_ssc = NULL, fixest_vcov = NULL,
                           poisson_model = NULL) {
  vce <- match.arg(vce); test <- match.arg(test); vcov_method <- match.arg(vcov_method)
  for (flag in list(eform, show_auxiliary, auxiliary_tests, show_deviance, only_test, small, diagnostics))
    if (!is.logical(flag) || length(flag) != 1 || is.na(flag)) stop("Display flags must be TRUE or FALSE.", call. = FALSE)
  if (only_test && test == "none") stop("only_test=TRUE requires a global test; test='none' disables it.", call. = FALSE)
  if (show_deviance && !inherits(model, "glm"))
    stop("show_deviance is available only for supported glm models.", call. = FALSE)
  iv <- d$type %in% c("iv_2sls", "fixest_iv")
  is_fixest <- inherits(model,"fixest")
  if (small && !iv) stop("small is for IV/2SLS models only.",call.=FALSE)
  if ((!is.null(fixest_ssc) || !is.null(fixest_vcov)) && !is_fixest)
    stop("fixest_ssc/fixest_vcov require a fixest model.",call.=FALSE)
  if (!is.null(vcov) && (!is.null(fixest_ssc) || !is.null(fixest_vcov)))
    stop("Custom vcov cannot be combined with fixest covariance controls.",call.=FALSE)
  if (diagnostics && (d$type != "iv_2sls" || vce != "standard" || !is.null(vcov) || !is.null(cluster)))
    stop("diagnostics currently reports conventional ivreg diagnostics only; not robust/cluster or fixest diagnostics.",call.=FALSE)
  if (iv && !is.null(null_model)) stop("null_model/LR is not defined for this 2SLS summary; use test_terms for Wald restrictions.",call.=FALSE)
  if (eform && !d$eform) stop("Exponentiation is not supported for this model/link.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1 || !is.finite(level)) stop("Invalid confidence level.", call. = FALSE)
  if (level > 1) level <- level/100
  if (level <= 0 || level >= 1) stop("level must be between 0 and 1.", call. = FALSE)
  if (!is.numeric(digits) || length(digits) != 1 || !is.finite(digits) || digits < 0 ||
      digits > 12 || digits != floor(digits)) stop("digits must be an integer from 0 to 12.", call. = FALSE)
  for (entry in list(correction, bic_n)) if (!is.null(entry) &&
      (!is.numeric(entry) || length(entry) != 1 || !is.finite(entry) || entry <= 0))
    stop("correction and bic_n must be positive finite numbers.", call. = FALSE)
  if (!is.null(df) && (!is.numeric(df) || length(df) != 1 || is.na(df) || df <= 0))
    stop("df must be positive or Inf.", call. = FALSE)
  if (!is.list(hessian_control) || (length(hessian_control) && is.null(names(hessian_control))))
    stop("hessian_control must be a named list.", call. = FALSE)
  if (!is.null(vcov_label) && (!is.character(vcov_label) || length(vcov_label) != 1 || is.na(vcov_label)))
    stop("vcov_label must be one character string.", call. = FALSE)
  if (!is.null(vcov) && (vce != "standard" || !is.null(hc_type) || !is.null(correction) ||
                         !is.null(cluster) || vcov_method != "auto"))
    stop("Supply vcov alone, without covariance-estimator options.", call. = FALSE)
  if (!is.null(hc_type) && (!d$type %in% c("linear", "logit", "probit", "poisson", "fixest_linear") || vce != "robust"))
    stop("hc_type requires vce='robust' and lm/binomial/Poisson glm or fixest OLS.", call. = FALSE)
  if (vce == "standard" && (!is.null(cluster) || !is.null(correction)))
    stop("cluster/correction require robust or cluster inference.", call. = FALSE)
  if (length(hessian_control) && d$type != "nested") stop("hessian_control is for nested logit only.", call. = FALSE)
  if (vcov_method == "opg" && vce != "standard") stop("opg is a separate estimator; do not combine it with robust/cluster.", call. = FALSE)
  linear <- d$type == "linear"
  b <- d$b; keep <- !is.na(b); nm <- names(b)[keep]
  if (!length(nm) || any(!is.finite(b[keep]))) stop("Finite estimable coefficients are required.", call. = FALSE)
  if ((linear || iv || d$type == "fixest_linear") && d$residual_df <= 0) stop("Positive residual degrees of freedom are required.", call. = FALSE)
  notes <- d$notes; altered <- !is.null(vcov) || vce != "standard" || !is.null(model$naive.var)
  method <- if (vcov_method == "auto") {
    if (d$type %in% c("ordinal_logit", "ordinal_probit", "nested", "negbin")) "observed" else "model"
  } else vcov_method
  detail <- NULL
  extension <- NULL
  if (isTRUE(d$extension)) {
    extension <- .ext_extension_covariance(model,d,vce,vcov,hc_type,correction,cluster,small,df,
                                            fixest_ssc,fixest_vcov,vcov_method)
    V <- extension$vcov; label <- extension$label; altered <- extension$altered
  }
  if (method == "observed" && is.null(vcov)) {
    if (d$type %in% c("ordinal_logit", "ordinal_probit")) detail <- .ext_polr_information(model)
    else if (d$type == "nested") detail <- .ext_nested_information(model, hessian_control)
    else if (d$type == "negbin") detail <- .ext_negbin_information(model)
    else stop("observed is available for polr and nested logit; other models use their native covariance.", call. = FALSE)
  }
  if (!isTRUE(d$extension)) {
  if (!is.null(vcov)) { V <- vcov; label <- "Custom" } else {
    if (d$type == "negbin" && is.null(detail))
      stop("NB2 needs joint observed information; vcov_method='model' does not include dispersion. Use auto/observed or explicit vcov.",call.=FALSE)
    V <- if (!is.null(detail)) detail$vcov else stats::vcov(model)
    label <- if (altered && !is.null(model$naive.var)) "Model (robust covariance)" else "Standard"
    if (method == "opg") {
      S <- .ext_scores(model, d); V <- solve(crossprod(S)); label <- "OPG"; altered <- TRUE
    }
  }
  G <- NULL
  if (vce != "standard") {
    if (d$type %in% c("logit", "probit") && vce == "robust" && is.null(hc_type)) {
      y <- stats::model.response(d$sample)
      binary <- if (is.factor(y)) nlevels(y) == 2 else !is.matrix(y) && all(y %in% c(0, 1))
      if (!binary || any(model$prior.weights != 1)) stop("Automatic binomial robust correction requires individual unweighted binary outcomes; supply vcov or explicit hc_type.", call. = FALSE)
    }
    if (is.null(cluster) && d$type == "conditional") cluster <- d$groups
    if (vce == "cluster" && is.null(cluster)) stop("Supply a grouping vector aligned to the estimation sample.", call. = FALSE)
    if (!is.null(cluster)) {
      if (!is.atomic(cluster) || anyNA(cluster) || length(cluster) != length(d$row_ids))
        stop("cluster must be one vector aligned to fitted observations/choice cases, with no missing values.", call. = FALSE)
      if (!is.null(names(cluster))) {
        if (anyDuplicated(names(cluster)) || !setequal(names(cluster), d$row_ids))
          stop("Named cluster vector must match estimation row/case identifiers exactly.", call. = FALSE)
        cluster <- cluster[d$row_ids]
      }
      if (!is.null(hc_type)) stop("Use hc_type for independent robust errors; clustering uses case scores.", call. = FALSE)
      S <- if (!is.null(detail)) detail$scores else .ext_scores(model, d)
      if (nrow(S) != length(cluster)) stop("Score rows do not match cluster vector.", call. = FALSE)
      G <- length(unique(cluster)); if (G < 2) stop("At least two clusters are required.", call. = FALSE)
      factor <- if (!is.null(correction)) correction else G/(G-1)*
        if (linear) (d$n-1)/d$residual_df else 1
      meat <- crossprod(rowsum(S, cluster, reorder = FALSE))
      if (linear) {
        if (!requireNamespace("sandwich", quietly = TRUE)) stop("Clustered lm requires sandwich.", call. = FALSE)
        bread <- sandwich::bread(model)/d$n
      } else bread <- if (!is.null(model$naive.var)) {
        raw <- model$naive.var; dimnames(raw) <- list(names(b), names(b)); raw
      } else V
      V <- bread %*% meat %*% bread * factor
      label <- paste0("Cluster (G=", G, ")")
    } else if (d$type %in% c("linear", "logit", "probit", "poisson")) {
      if (!requireNamespace("sandwich", quietly = TRUE)) stop("Robust inference requires sandwich.", call. = FALSE)
      h <- if (is.null(hc_type)) if (linear) "HC1" else "HC0" else
        match.arg(hc_type, c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5"))
      factor <- if (!is.null(correction)) correction else if (!linear && is.null(hc_type)) d$n/(d$n-1) else 1
      V <- sandwich::vcovHC(model, type = h)*factor
      label <- if (!linear && is.null(hc_type) && is.null(correction)) "Robust (HC0 * N/(N-1))" else paste0("Robust (", h, ")")
    } else {
      S <- if (!is.null(detail)) detail$scores else .ext_scores(model, d)
      if (!is.null(d$weights) && any(d$weights != 1) && is.null(correction))
        stop("For weighted discrete robust inference, specify correction or supply vcov.", call. = FALSE)
      N <- nrow(S); if (N <= 1) stop("Robust inference requires multiple independent cases.", call. = FALSE)
      factor <- if (is.null(correction)) N/(N-1) else correction
      V <- V %*% crossprod(S) %*% V * factor
      label <- "Robust (case scores)"
    }
  }
  } else G <- extension$clusters
  if (!is.null(vcov_label)) label <- vcov_label
  V <- .ext_validate_vcov(V, nm)
  infer_df <- if (!is.null(extension)) extension$df else if (!is.null(df)) df else if (linear) if (is.null(G)) d$residual_df else G-1 else Inf
  se <- rep(NA_real_, length(b)); se[keep] <- sqrt(diag(V))
  statistic <- b/se; alpha <- (1-level)/2
  critical <- if (is.finite(infer_df)) stats::qt(1-alpha, infer_df) else stats::qnorm(1-alpha)
  p <- if (is.finite(infer_df)) 2*stats::pt(-abs(statistic), infer_df) else 2*stats::pnorm(-abs(statistic))
  aux <- names(b) %in% d$auxiliary
  if (!auxiliary_tests) { statistic[aux] <- NA_real_; p[aux] <- NA_real_ }
  transform <- eform & !aux
  lower <- b-critical*se; upper <- b+critical*se
  estimate <- b; displayed_se <- se
  estimate[transform] <- exp(b[transform]); displayed_se[transform] <- estimate[transform]*se[transform]
  lower[transform] <- exp(lower[transform]); upper[transform] <- exp(upper[transform])
  tab <- data.frame(estimate = estimate, std.error = displayed_se, statistic = statistic, p.value = p,
                    conf.low = lower, conf.high = upper, omitted = !keep,
                    auxiliary = aux, group = d$group, term = d$term,
                    row.names = names(b), check.names = FALSE)
  k <- d$k; n <- if (is.null(bic_n)) d$n else bic_n
  fit <- list(n = d$n, parameters = k, log_likelihood = d$ll,
              AIC = -2*d$ll+2*k, BIC = -2*d$ll+log(n)*k, bic_n = n,
              null_log_likelihood = d$ll0)
  if (inherits(model, "glm")) {
    fit$null_deviance <- model$null.deviance
    fit$null_df <- model$df.null
    fit$residual_deviance <- model$deviance
    fit$residual_df <- model$df.residual
  }
  if (!is.null(d$fit_extra)) fit <- c(fit, d$fit_extra)
  if (!is.null(null_model)) {
    nd <- if (inherits(null_model,"negbin")) .ext_negbin(null_model) else
      if (inherits(null_model,"fixest")) .ext_fixest(null_model) else
      if (inherits(null_model, "glm")) .ext_glm(null_model) else
      if (inherits(null_model, "lm")) .ext_lm(null_model) else
      if (inherits(null_model, "multinom")) .ext_multinom(null_model) else
      if (inherits(null_model, "polr")) .ext_polr(null_model) else
      if (inherits(null_model, "clogit")) .ext_clogit(null_model) else
      if (inherits(null_model, "mlogit")) .ext_mlogit(null_model) else
      stop("Unsupported null_model.", call. = FALSE)
    if (nd$type != d$type || !identical(nd$row_ids, d$row_ids) || nd$n != d$n)
      stop("null_model must use the same model family and identical estimation sample/order.", call. = FALSE)
    if (d$type == "conditional" && model$method != null_model$method)
      stop("null_model must use the same conditional likelihood method.", call. = FALSE)
    if (d$type == "nested" && !identical(model$nests, null_model$nests))
      stop("null_model must use the same nesting structure.", call. = FALSE)
    full_response <- if (is_fixest) stats::model.matrix(model,type="lhs") else stats::model.response(stats::model.frame(model))
    null_response <- if (is_fixest) stats::model.matrix(null_model,type="lhs") else stats::model.response(stats::model.frame(null_model))
    if (!isTRUE(all.equal(full_response, null_response))) stop("null_model must have the same response.", call. = FALSE)
    full_w <- if (is_fixest) model$weights else stats::model.weights(stats::model.frame(model))
    null_w <- if (is_fixest) null_model$weights else stats::model.weights(stats::model.frame(null_model))
    if (!isTRUE(all.equal(full_w, null_w))) stop("null_model must use the same weights.", call. = FALSE)
    if (is_fixest && !identical(model$fixef_id,null_model$fixef_id))
      stop("fixest null_model must retain identical absorbed effects; LR tests only the additional slopes.",call.=FALSE)
    if (is_fixest && !isTRUE(all.equal(model$offset,null_model$offset)))
      stop("fixest null_model must retain identical offsets/exposure.",call.=FALSE)
    if (d$type %in% c("poisson","negbin")) {
      fo <- stats::model.offset(stats::model.frame(model)); no <- stats::model.offset(stats::model.frame(null_model))
      if (!isTRUE(all.equal(fo,no))) stop("Count null_model must retain the same offsets/exposure.",call.=FALSE)
    }
    fit$null_log_likelihood <- nd$ll; d$null_df <- nd$k
    notes <- c(notes, "Null statistics use the supplied null_model; LR validity requires a nested specification.")
  }
  fit$pseudo_r_squared <- if (is.finite(fit$null_log_likelihood) && fit$null_log_likelihood != 0)
    1-d$ll/fit$null_log_likelihood else NA_real_
  joint <- if (is.null(test_terms)) d$global else {
    if (!is.character(test_terms) || anyDuplicated(test_terms) || !all(test_terms %in% nm))
      stop("test_terms must contain distinct exact names of estimable parameters.", call. = FALSE)
    test_terms
  }
  joint <- intersect(joint, nm)
  if (!is.null(test_terms) && test == "LR") stop("test_terms defines a Wald test; use test='Wald'.", call. = FALSE)
  if (test == "auto") test <- if (linear || altered || !is.null(test_terms) ||
      isTRUE(d$force_wald) || d$type %in% c("alternative_specific", "nested")) "Wald" else "LR"
  if (test == "LR" && (linear || iv || d$type == "fixest_linear" || altered)) stop("LR is not available with robust/custom covariance, OLS or 2SLS; use Wald.", call. = FALSE)
  if (test == "LR" && is_fixest && is.null(null_model)) stop("fixest LR requires an explicit compatible null_model with the same fixed effects.",call.=FALSE)
  global <- list(type = "none", statistic = NA_real_, df = 0L, denominator_df = NULL, p.value = NA_real_, terms = joint)
  if (test == "Wald" && length(joint)) {
    stat <- tryCatch(as.numeric(crossprod(b[joint], solve(V[joint, joint, drop = FALSE], b[joint]))),
                     error = function(e) NA_real_)
    global$type <- if (is.finite(infer_df)) "F" else "Wald chi2"
    global$df <- length(joint); global$denominator_df <- if (is.finite(infer_df)) infer_df else NULL
    global$statistic <- if (is.finite(infer_df)) stat/length(joint) else stat
    global$p.value <- if (is.finite(infer_df)) stats::pf(global$statistic, length(joint), infer_df, lower.tail = FALSE) else
      stats::pchisq(stat, length(joint), lower.tail = FALSE)
    if (is.na(stat)) notes <- c(notes, "Joint test unavailable: covariance matrix is singular.")
  }
  if (test == "LR") {
    if (!is.finite(fit$null_log_likelihood)) stop("Null likelihood is unavailable for this specification; supply null_model or use test='Wald'.", call. = FALSE)
    global$type <- "LR chi2"; global$df <- k-d$null_df
    LR <- 2*(d$ll-fit$null_log_likelihood)
    if (LR < -1e-6 || global$df < 0) stop("Invalid LR comparison; check nesting, sample and optimization.", call. = FALSE)
    if (global$df > 0) {
      global$statistic <- max(0, LR); global$p.value <- stats::pchisq(global$statistic, global$df, lower.tail = FALSE)
    }
  }
  iia <- NULL
  if (!is.null(iia_model)) {
    if (d$type != "nested" || altered || !inherits(iia_model, "mlogit") || !is.null(iia_model$nests))
      stop("iia_model requires conventional nested inference and an ordinary mlogit.", call. = FALSE)
    restricted <- .ext_mlogit(iia_model)
    if (!identical(d$ix, restricted$ix) ||
        !isTRUE(all.equal(stats::model.matrix(model), stats::model.matrix(iia_model), check.attributes = FALSE)) ||
        !isTRUE(all.equal(as.logical(stats::model.response(model$model)),
                         as.logical(stats::model.response(iia_model$model)))))
      stop("iia_model must have identical alternatives, cases, choices and utility regressors.", call. = FALSE)
    if (!is.null(stats::model.weights(model$model)) || !is.null(stats::model.weights(iia_model$model)))
      stop("The nested IIA LR comparison currently requires unweighted models.", call. = FALSE)
    difference <- d$k-restricted$k
    LR <- 2*(d$ll-restricted$ll)
    if (difference != length(d$auxiliary) || LR < -1e-6)
      stop("Invalid nested-versus-MNL comparison; check specifications and convergence.", call. = FALSE)
    iia <- list(type = "LR chi2", statistic = max(0, LR), df = difference,
                p.value = stats::pchisq(max(0, LR), difference, lower.tail = FALSE))
  }
  if (linear) {
    native <- summary(model); fit$residual_df <- d$residual_df
    fit$r_squared <- native$r.squared; fit$adj_r_squared <- native$adj.r.squared; fit$root_mse <- native$sigma
    if (!altered && is.null(model$offset)) {
      f <- model$fitted.values; w <- model$weights
      if (is.null(w)) w <- rep(1, length(f))
      intercept <- attr(model$terms, "intercept"); center <- if (intercept) stats::weighted.mean(f, w) else 0
      SS <- c(sum(w*(f-center)^2), sum(w*model$residuals^2)); SS <- c(SS, sum(SS))
      dfs <- c(model$rank-intercept, d$residual_df, model$rank-intercept+d$residual_df)
      fit$source <- data.frame(SS = SS, df = dfs, MS = ifelse(dfs > 0, SS/dfs, NA_real_),
                              row.names = c("Model", "Residual", "Total"))
    }
    if (attr(model$terms, "intercept") == 0) notes <- c(notes, "R-squared is uncentered because the model has no intercept.")
  }
  if (iv || d$type == "fixest_linear") {
    fit$root_mse <- sqrt(d$rss/if (iv && !small) d$n else d$residual_df)
    if (iv) fit$AIC <- fit$BIC <- fit$log_likelihood <- fit$null_log_likelihood <- fit$pseudo_r_squared <- NA_real_
  }
  if (diagnostics) {
    if (!requireNamespace("ivreg",quietly=TRUE)) stop("ivreg package is required for diagnostics.",call.=FALSE)
    fit$iv_diagnostics <- summary(model,diagnostics=TRUE)$diagnostics
    notes <- c(notes,"IV diagnostic table uses conventional ivreg first-stage F, Wu-Hausman and Sargan tests; not robust weak-IV inference.")
  }
  if (d$type %in% c("negbin","fixest_negbin")) {
    pos <- match("lnalpha",rownames(tab)); fit$alpha <- exp(b[pos])
    fit$alpha_se <- fit$alpha*se[pos]
    fit$alpha_ci <- exp(c(b[pos]-critical*se[pos],b[pos]+critical*se[pos]))
  }
  dispersion_test <- NULL
  if (!is.null(poisson_model)) {
    if (!d$type %in% c("negbin","fixest_negbin") || altered)
      stop("poisson_model requires an NB2 model with conventional covariance.",call.=FALSE)
    if (is_fixest) {
      pd <- .ext_fixest(poisson_model)
      compatible <- pd$type == "fixest_poisson" && identical(model$fixef_id,poisson_model$fixef_id) &&
        isTRUE(all.equal(stats::model.matrix(model,type="lhs"),stats::model.matrix(poisson_model,type="lhs"))) &&
        isTRUE(all.equal(stats::model.matrix(model),stats::model.matrix(poisson_model))) &&
        isTRUE(all.equal(model$offset,poisson_model$offset))
    } else {
      if (!inherits(poisson_model,"glm") || inherits(poisson_model,"negbin")) stop("Supply a fitted Poisson glm.",call.=FALSE)
      pd <- .ext_glm(poisson_model)
      compatible <- pd$type == "poisson" &&
        isTRUE(all.equal(stats::model.response(d$sample),stats::model.response(pd$sample))) &&
        isTRUE(all.equal(stats::model.matrix(model),stats::model.matrix(poisson_model))) &&
        isTRUE(all.equal(stats::model.offset(d$sample),stats::model.offset(pd$sample)))
    }
    if (!compatible || !identical(d$row_ids,pd$row_ids)) stop("Poisson comparison must use identical sample, regressors, offsets and fixed effects.",call.=FALSE)
    LR <- 2*(d$ll-pd$ll)
    if (LR < -1e-6) stop("NB2 has lower likelihood than Poisson; inspect convergence/boundary.",call.=FALSE)
    LR <- max(0,LR)
    dispersion_test <- list(type="LR chibar2(01)",statistic=LR,df=1,
                             p.value=if (LR == 0) 1 else .5*stats::pchisq(LR,1,lower.tail=FALSE))
  }
  if (eform) notes <- c(notes, if (d$type %in% c("multinomial", "multinomial_choice"))
    "Exponentiated slopes are relative risk ratios versus the base outcome." else
    if (d$type %in% c("poisson","negbin","fixest_poisson","fixest_negbin"))
      "Exponentiated slopes are incidence-rate ratios (IRR); dispersion stays on its original scale." else
      "Exponentiated slopes are odds ratios; auxiliary parameters remain on their original scale.")
  if (any(!keep)) notes <- c(notes, paste("Omitted (aliased):", paste(names(b)[!keep], collapse = ", ")))
  if (any(!is.finite(tab$estimate[keep]) | !is.finite(tab$std.error[keep]) | !is.finite(tab$conf.high[keep])))
    notes <- c(notes, "Some displayed estimates or intervals overflow; inspect the original coefficient scale.")
  if (!is.finite(fit$null_log_likelihood) && !linear && !iv && d$type != "fixest_linear") notes <- c(notes, "Null likelihood and pseudo-R2 unavailable; supply a compatible null_model.")
  result <- structure(list(model_type = d$type, model = model, coefficients = tab, vcov = V,
                           vce = label, global_test = global, fit = fit, eform = eform, level = level,
                           digits = as.integer(digits), notes = unique(notes),
                           base = d$base, inference_df = infer_df, show_auxiliary = show_auxiliary,
                           vcov_method = method, clusters = G, iia_test = iia,
                           show_deviance = show_deviance, only_test = only_test,
                           small = small, covariance_details = if (!is.null(extension)) extension$metadata else NULL,
                           dispersion_test = dispersion_test), class = "ext_summary")
  print(result)
  invisible(result)
}


#' @rdname ext_summary
#' @export
print.ext_summary <- function(x, digits = x$digits, only_test = x$only_test, ...) {
  if (!is.null(only_test) && (!is.logical(only_test) || length(only_test) != 1 || is.na(only_test)))
    stop("only_test must be TRUE or FALSE.", call. = FALSE)
  title <- switch(x$model_type, linear = "Linear regression", logit = "Logistic regression",
                  iv_2sls = "Instrumental-variables 2SLS regression",
                  poisson = "Poisson regression", negbin = "Negative binomial regression (NB2)",
                  fixest_linear = "Linear regression with fixest",
                  fixest_iv = "Instrumental-variables 2SLS with fixest",
                  fixest_poisson = "Poisson regression with fixest",
                  fixest_negbin = "Negative binomial regression with fixest (NB2)",
                  fixest_logit = "Logistic regression with fixest", fixest_probit = "Probit regression with fixest",
                  probit = "Probit regression", multinomial = "Multinomial logistic regression",
                  multinomial_choice = "Multinomial logistic regression",
                  ordinal_logit = "Ordered logistic regression", ordinal_probit = "Ordered probit regression",
                  conditional = "Conditional logistic regression",
                  alternative_specific = "Alternative-specific conditional logit", nested = "Nested logit regression")
  fmt <- function(value, places = digits) {
    if (!is.finite(value)) return(".")
    if ((value != 0 && abs(value) < 10^(-places)) || abs(value) >= 1e8)
      formatC(value, format = "e", digits = 2) else formatC(value, format = "f", digits = places)
  }
  field <- function(label, value) cat(sprintf("%-26s = %s\n", label, value))
  if (isTRUE(only_test)) {
    g <- x$global_test
    test_name <- switch(g$type,
                        "LR chi2" = "Likelihood-ratio test (LR)",
                        "Wald chi2" = "Joint Wald test",
                        "F" = "Joint F test",
                        "Global test unavailable")
    cat("\n", test_name, "\n", sep = "")
    cat("Model: ", title, "\n", sep = "")
    if (g$df > 0 && g$type != "none") {
      label <- if (g$type == "F") sprintf("F(%d, %s)", g$df, format(g$denominator_df, trim = TRUE))
               else sprintf("%s(%d)", g$type, g$df)
      field(label, fmt(g$statistic, 2))
      field(if (g$type == "F") "Prob > F" else "Prob > chi2", fmt(g$p.value, 4))
      if (!is.finite(g$statistic)) cat("Note: Joint test unavailable; inspect the fitted model and covariance matrix.\n")
    } else cat("Global test unavailable: no tested degrees of freedom or test disabled.\n")
    cat("\n")
    return(invisible(x))
  }
  cat("\nExtended Summary - ", title, "\n", sep = "")
  formula <- if (x$model_type == "conditional" && !is.null(x$model$userCall))
    x$model$userCall[[2L]] else stats::formula(x$model)
  cat("Formula: ", paste(deparse(formula), collapse = " "), "\n", sep = "")
  cat("\n")
  header <- c(sprintf("%-26s = %s", "Number of obs", as.character(x$fit$n)),
              sprintf("%-26s = %s", "VCE", x$vce))
  add_field <- function(label, value) {
    header <<- c(header, sprintf("%-26s = %s", label, value))
  }
  if (!is.null(x$fit$cases)) {
    add_field("Number of cases", as.character(x$fit$cases))
    if (!is.null(x$fit$informative_cases)) add_field("Informative cases", as.character(x$fit$informative_cases))
    add_field("Alternatives/case (min)", fmt(x$fit$alternatives[1], 0))
    add_field("Alternatives/case (avg)", fmt(x$fit$alternatives[2], 1))
    add_field("Alternatives/case (max)", fmt(x$fit$alternatives[3], 0))
  }
  if (length(x$base)) add_field("Base outcome/alternative", paste(x$base, collapse = ", "))
  if (length(x$fit$fixed_effects)) add_field("Fixed effects (levels)",paste(paste0(names(x$fit$fixed_effects),": ",x$fit$fixed_effects),collapse=", "))
  g <- x$global_test
  if (g$df > 0) {
    test_label <- if (g$type == "F") sprintf("F(%d, %s)", g$df, format(g$denominator_df, trim = TRUE))
                  else sprintf("%s(%d)", g$type, g$df)
    add_field(test_label, fmt(g$statistic, 2))
    add_field(if (g$type == "F") "Prob > F" else "Prob > chi2", fmt(g$p.value, 4))
  }
  if (x$model_type %in% c("linear","iv_2sls","fixest_linear","fixest_iv")) {
    add_field("R-squared", fmt(x$fit$r_squared, 4))
    add_field("Adjusted R-squared", fmt(x$fit$adj_r_squared, 4))
    add_field("Root MSE", fmt(x$fit$root_mse))
    if (!is.null(x$fit$within_r_squared) && is.finite(x$fit$within_r_squared))
      add_field("Within R-squared",fmt(x$fit$within_r_squared,4))
  } else {
    add_field("Log likelihood (Model)", fmt(x$fit$log_likelihood))
    add_field("Log likelihood (Null)", fmt(x$fit$null_log_likelihood))
    add_field("Pseudo R2 (McFadden)", fmt(x$fit$pseudo_r_squared, 4))
  }
  if (!is.null(x$fit$source)) {
    source <- x$fit$source
    ss <- vapply(source$SS, fmt, character(1))
    ms <- vapply(source$MS, fmt, character(1))
    ss_width <- max(14L, nchar(ss))
    ms_width <- max(14L, nchar(ms))
    source_line <- function(label, ss, df, ms) {
      sprintf("%-10s | %*s %6s %*s", label, ss_width, ss, df, ms_width, ms)
    }
    left <- c(source_line("Source", "SS", "df", "MS"),
              strrep("-", 10L + 3L + ss_width + 7L + 1L + ms_width))
    left <- c(left, source_line("Model", ss[1], source$df[1], ms[1]),
              source_line("Residual", ss[2], source$df[2], ms[2]), left[2],
              source_line("Total", ss[3], source$df[3], ms[3]))
    width <- max(nchar(left))
    for (i in seq_len(max(length(left), length(header)))) {
      cat(sprintf("%-*s    %s\n", width,
                  if (i <= length(left)) left[i] else "",
                  if (i <= length(header)) header[i] else ""))
    }
  } else {
    cat(paste(header, collapse = "\n"), "\n", sep = "")
  }
  tab <- x$coefficients
  if (!x$show_auxiliary) tab <- tab[!tab$auxiliary, , drop = FALSE]
  if (x$model_type %in% c("alternative_specific", "nested", "multinomial_choice")) {
    groups <- c("Generic slopes", setdiff(unique(tab$group), c("Generic slopes", "Dissimilarity")), "Dissimilarity")
    tab <- tab[order(match(tab$group, groups), tab$term == "(Intercept)"), , drop = FALSE]
  }
  if (is.null(tab$term)) tab$term <- rownames(tab)
  name_width <- max(12L, nchar(tab$term))
  number_width <- max(12L, digits + 6L)
  separator <- strrep("-", name_width + 6L * (number_width + 1L) + 9L)
  cat("\n", separator, "\n", sep = "")
  labels <- c(if (x$eform) if (x$model_type %in% c("multinomial", "multinomial_choice")) "RRR" else
              if (x$model_type %in% c("poisson","negbin","fixest_poisson","fixest_negbin")) "IRR" else "Odds Ratio" else "Coefficient", "Std. Error",
              if (is.finite(x$inference_df)) "t" else "z",
              if (is.finite(x$inference_df)) "P>|t|" else "P>|z|", "CI Lower", "CI Upper")
  cat(sprintf("%-*s |", name_width, "Variable"))
  header_values <- sprintf("%*s", number_width, labels)
  cat(paste(c(header_values[1:4], " Sig.", header_values[5:6]), collapse = " "), "\n", sep = "")
  cat(separator, "\n", sep = "")
  for (i in seq_len(nrow(tab))) {
    if (nzchar(tab$group[i]) && (i == 1L || tab$group[i] != tab$group[i-1L])) {
      cat(tab$group[i], "\n", sep = "")
    }
    values <- c(fmt(tab$estimate[i]), fmt(tab$std.error[i]), fmt(tab$statistic[i], 2),
                fmt(tab$p.value[i], 3), fmt(tab$conf.low[i]), fmt(tab$conf.high[i]))
    cat(sprintf("%-*s |", name_width, tab$term[i]))
    significance <- if (is.na(tab$p.value[i])) "" else
      as.character(stats::symnum(tab$p.value[i], corr = FALSE, na = FALSE,
                   cutpoints = c(0, .001, .01, .05, .1, 1),
                   symbols = c("***", "**", "*", ".", " ")))
    formatted <- sprintf("%*s", number_width, values)
    cat(paste(c(formatted[1:4], sprintf("%5s", significance), formatted[5:6]),
              collapse = " "), "\n", sep = "")
  }
  cat(separator, "\n", sep = "")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  if (isTRUE(x$show_deviance)) {
    cat(sprintf("\n    Null deviance: %s on %s degrees of freedom\n",
                fmt(x$fit$null_deviance, 3), x$fit$null_df))
    cat(sprintf("Residual deviance: %s on %s degrees of freedom\n\n",
                fmt(x$fit$residual_deviance, 3), x$fit$residual_df))
  }
  field("Confidence level", paste0(fmt(100 * x$level, 1), "%"))
  if (!x$model_type %in% c("iv_2sls","fixest_iv")) {
    field("AIC", fmt(x$fit$AIC))
    field("BIC", fmt(x$fit$BIC))
  }
  if (!is.null(x$fit$alpha) && x$show_auxiliary) {
    field("alpha (NB2 dispersion)",fmt(x$fit$alpha))
    field("Std. Error (alpha)",fmt(x$fit$alpha_se))
    field("CI (alpha)",paste(vapply(x$fit$alpha_ci,fmt,character(1)),collapse="  "))
  }
  if (!is.null(x$dispersion_test)) {
    field("LR alpha=0: chibar2(01)",fmt(x$dispersion_test$statistic,2))
    field("Prob >= chibar2",fmt(x$dispersion_test$p.value,4))
  }
  if (!is.null(x$fit$iv_diagnostics)) {
    cat("\nIV diagnostics (conventional):\n")
    print(x$fit$iv_diagnostics,digits=digits)
  }
  if (x$model_type %in% c("conditional", "alternative_specific", "nested") || x$fit$bic_n != x$fit$n)
    field("N used for BIC", as.character(x$fit$bic_n))
  if (!is.null(x$iia_test)) {
    field(sprintf("LR IIA: tau=1 (%d df)", x$iia_test$df), fmt(x$iia_test$statistic, 2))
    field("Prob > chi2 (IIA)", fmt(x$iia_test$p.value, 4))
  }
  for (note in x$notes) cat("Note: ", note, "\n", sep = "")
  cat("\n")
  invisible(x)
}
