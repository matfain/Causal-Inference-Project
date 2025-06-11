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

bootstrap_iptw_ate <- function(df,
                               confounders,
                               treatment_col,
                               outcome_col,
                               lasso = FALSE,
                               n_bootstrap = 200,
                               stabilized = TRUE,
                               clip = TRUE,
                               clip_quantiles = c(0.05, 0.95),
                               plot_title = "Bootstrap Distribution of IPTW‐ATE") {
  # parallel setup
  library(doParallel)
  library(foreach)
  cores <- parallel::detectCores() - 1
  cl    <- makeCluster(cores)
  registerDoParallel(cl)
  
  # run bootstrap in parallel, exporting the model & ATE functions
  ates <- foreach(
    b        = seq_len(n_bootstrap),
    .combine = c,
    .packages = c("glmnet"),
    .export   = c("model_ps_score", "extract_ps", "iptw_ate", "iptw_weights")
  ) %dopar% {
    # 1) resample
    idx  <- sample(seq_len(nrow(df)), replace = TRUE)
    df_b <- df[idx, , drop = FALSE]
    
    # 2) fit PS
    ps_fit <- model_ps_score(df_b, confounders, treatment_col, lasso)
    ps_vec <- extract_ps(ps_fit, df_b, confounders, treatment_col)
    
    # 3) compute ATE
    iptw_ate(
      df_b,
      outcome_col,
      treatment_col,
      ps_vec,
      stabilized     = stabilized,
      clip            = clip,
      clip_quantiles  = clip_quantiles
    )$ATE
  }
  
  stopCluster(cl)
  
  # percentile CI & SD
  ci     <- quantile(ates, probs = c(0.025, 0.975), na.rm = TRUE)
  sd_est <- sd(ates, na.rm = TRUE)
  
  # print results
  cat(sprintf("95%% percentile CI: [%.3f, %.3f]\n", ci[1], ci[2]))
  cat(sprintf("SD of ATE estimator: %.3f\n", sd_est))
  
  # histogram with CI lines
  library(ggplot2)
  df_plot <- data.frame(ATE = ates)
  p <- ggplot(df_plot, aes(x = ATE)) +
    geom_histogram(bins = 30, color = "black", fill = "plum") +
    geom_vline(xintercept = ci, linetype = "dashed", size = 1) +
    labs(
      title = plot_title,
      x     = "ATE",
      y     = "Frequency"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  print(p)
  
  invisible(list(ates = ates, ci = ci, sd = sd_est))
}
