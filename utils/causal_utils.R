# function to accept a df, confounders, treatment, boolean value for lasso regularization.
# Return logistic regression model with lasso regularization if boolean is true, otherwise return a model without lasso regularization.
# If lasso is true then use CV to find the best lambda.

model_ps_score <- function(df, confounders, treatment, lasso = FALSE) {
  fmla <- as.formula(paste(treatment, "~", paste(confounders, collapse = " + ")))
  
  if (!lasso) {
    # plain logistic‐glm
    return(glm(fmla, data = df, family = binomial))
  }
  
  # build design matrix (drop intercept col)
  X <- model.matrix(fmla, data = df)[, -1, drop = FALSE]
  y <- df[[treatment]]
  if (is.factor(y)) y <- as.numeric(y) - 1
  
  cvfit <- glmnet::cv.glmnet(x = X, y = y, family = "binomial", alpha = 1)
  lam0  <- cvfit$lambda.min
  
  # refit on full data at that single lambda
  glmnet::glmnet(
    x      = X,
    y      = y,
    family = "binomial",
    alpha  = 1,
    lambda = lam0
  )
}

extract_ps <- function(fit, df = NULL, confounders = NULL, treatment = NULL) {
  if (inherits(fit, "glm")) {
    # fitted values are the propensity scores
    return(fitted(fit))
  }
  
  if (inherits(fit, "glmnet")) {
    if (is.null(df) || is.null(confounders) || is.null(treatment)) {
      stop("For glmnet fits you must supply df, confounders, and treatment.")
    }
    # rebuild design matrix (drop intercept)
    fmla <- as.formula(paste(treatment, "~", paste(confounders, collapse = " + ")))
    X <- model.matrix(fmla, data = df)[, -1, drop = FALSE]
    # predict at the single λ stored in fit$lambda
    return(as.vector(predict(fit, newx = X, type = "response")))
  }
  
  stop("`fit` must be either a glm or a glmnet object.")
}

iptw_weights <- function(df,
                         ps_vec,
                         treatment_col,
                         stabilized     = FALSE,
                         clip            = FALSE,
                         clip_quantiles = c(0.05, 0.95)) {
  stopifnot(length(ps_vec) == nrow(df))
  
  # 1) extract & coerce treatment to 0/1
  trt <- df[[treatment_col]]
  if (is.factor(trt)) trt <- as.numeric(trt) - 1
  trt <- as.numeric(trt)
  
  # 2) propensity scores
  ps <- as.numeric(ps_vec)
  
  # 3) compute weights
  if (stabilized) {
    p_t <- mean(trt, na.rm = TRUE)
    w   <- ifelse(trt == 1,
                  p_t       / ps,
                  (1 - p_t) / (1 - ps))
  } else {
    w <- ifelse(trt == 1,
                1 / ps,
                1 / (1 - ps))
  }
  
  # 4) optional clipping
  if (clip) {
    q <- quantile(w, probs = clip_quantiles, na.rm = TRUE)
    w[w <  q[1]] <- q[1]
    w[w >  q[2]] <- q[2]
  }
  
  w
}

iptw_ate <- function(df,
                     outcome_col,
                     treatment_col,
                     ps_vec,
                     stabilized     = FALSE,
                     clip            = FALSE,
                     clip_quantiles = c(0.05, 0.95)) {
  # get weights (with clipping if requested)
  w <- iptw_weights(df,
                    ps_vec,
                    treatment_col,
                    stabilized     = stabilized,
                    clip            = clip,
                    clip_quantiles = clip_quantiles)
  
  # extract vars
  y   <- df[[outcome_col]]
  trt <- df[[treatment_col]]
  if (is.factor(trt)) trt <- as.numeric(trt) - 1
  trt <- as.numeric(trt)
  
  # weighted means
  mu1 <- sum(w * trt     * y, na.rm = TRUE) / sum(w * trt, na.rm = TRUE)
  mu0 <- sum(w * (1-trt) * y, na.rm = TRUE) / sum(w * (1-trt), na.rm = TRUE)
  ate <- mu1 - mu0
  
  list(
    weights      = w,
    mean_treated = mu1,
    mean_control = mu0,
    ATE          = ate
  )
}