# =======================================================================
#  Standardisation / G-Computation helpers
#  Place this file under utils/  and source() it in the notebook.
# =======================================================================

# -----------------------------------------------------------------------
# 1. Fit outcome model  P(Y = 1 | A, X)
# -----------------------------------------------------------------------
model_outcome <- function(df,
                          covariates,
                          treatment,
                          outcome,
                          model_type = c("logreg", "logreg_l1", "rf")) {
  model_type <- match.arg(model_type)
  
  # Build formula: outcome ~ treatment + covariates
  fmla <- as.formula(
    paste(outcome, "~", treatment, "+", paste(covariates, collapse = " + "))
  )
  
  # ----- Logistic GLM ---------------------------------------------------
  if (model_type == "logreg") {
    return(glm(fmla, data = df, family = binomial))
  }
  
  # ----- L1-penalised logistic via glmnet -------------------------------
  if (model_type == "logreg_l1") {
    X <- model.matrix(fmla, data = df)[, -1, drop = FALSE]   # drop intercept
    y <- df[[outcome]]
    cv <- glmnet::cv.glmnet(x = X, y = y,
                            family = "binomial",
                            alpha  = 1)                      # L1
    lam <- cv$lambda.min
    return(glmnet::glmnet(x = X, y = y,
                          family = "binomial",
                          alpha  = 1,
                          lambda = lam))
  }
  
  # ----- Random forest (ranger) ----------------------------------------
  if (model_type == "rf") {
    df_rf <- df %>%
      mutate(!!outcome := factor(.data[[outcome]], levels = c(0, 1)))
    
    return(ranger::ranger(
      formula       = fmla,
      data          = df_rf,
      probability   = TRUE,
      num.trees     = 500,
      mtry          = floor(sqrt(length(covariates) + 1)),   # simple default
      class.weights = as.numeric(1 / table(df_rf[[outcome]]))
    ))
  }
  
  stop("Unsupported model_type")
}

# -----------------------------------------------------------------------
# 2. Predict counterfactual risks  Ê[Y | A = a, X = x_i]
# -----------------------------------------------------------------------
predict_counterfactual <- function(fit,
                                   df,
                                   covariates,
                                   treatment,
                                   outcome) {
  
  # Helper to overwrite treatment assignment
  switch_trt <- function(dat, a_val) {
    dat2 <- dat
    dat2[[treatment]] <- a_val
    dat2
  }
  
  if (inherits(fit, "glm")) {
    p1 <- predict(fit, newdata = switch_trt(df, 1), type = "response")
    p0 <- predict(fit, newdata = switch_trt(df, 0), type = "response")
    
  } else if (inherits(fit, "glmnet")) {
    fmla <- as.formula(
      paste(outcome, "~", treatment, "+", paste(covariates, collapse = " + "))
    )
    X1 <- model.matrix(fmla, data = switch_trt(df, 1))[, -1, drop = FALSE]
    X0 <- model.matrix(fmla, data = switch_trt(df, 0))[, -1, drop = FALSE]
    p1 <- predict(fit, newx = X1, type = "response")[, 1]
    p0 <- predict(fit, newx = X0, type = "response")[, 1]
    
  } else if (inherits(fit, "ranger")) {
    p1 <- predict(fit, switch_trt(df, 1))$predictions[, "1"]
    p0 <- predict(fit, switch_trt(df, 0))$predictions[, "1"]
    
  } else {
    stop("Unsupported fit object")
  }
  
  tibble(y1_hat = p1, y0_hat = p0)
}

# -----------------------------------------------------------------------
# 3. Point-estimate of ATE via standardisation
# -----------------------------------------------------------------------
standardization_ate <- function(pred_df) {
  mu1 <- mean(pred_df$y1_hat)
  mu0 <- mean(pred_df$y0_hat)
  list(mu1 = mu1, mu0 = mu0, ATE = mu1 - mu0)
}

# -----------------------------------------------------------------------
# 4. Bootstrap for variance / CI
# -----------------------------------------------------------------------
bootstrap_standardization_ate <- function(df,
                                          covariates,
                                          treatment,
                                          outcome,
                                          model_type,
                                          n_bootstrap = 200,
                                          plot_title  = "Bootstrap ATE (standardisation)") {
  
  library(doParallel); library(foreach)
  
  cores <- parallel::detectCores() - 1
  cl    <- makeCluster(cores)
  registerDoParallel(cl)
  
  ates <- foreach(b = seq_len(n_bootstrap),
                  .combine  = c,
                  .packages = c("glmnet", "ranger", "tidyverse", "stats"),
                  .export   = c("model_outcome", "predict_counterfactual",
                                "standardization_ate")) %dopar% {
                                  
                                  idx  <- sample(seq_len(nrow(df)), replace = TRUE)
                                  df_b <- df[idx, , drop = FALSE]
                                  
                                  fit   <- model_outcome(df_b, covariates, treatment, outcome, model_type)
                                  preds <- predict_counterfactual(fit, df_b, covariates, treatment, outcome)
                                  standardization_ate(preds)$ATE
                                }
  
  stopCluster(cl)
  
  ci     <- quantile(ates, probs = c(0.025, 0.975))
  sd_est <- sd(ates)
  
  # Quick histogram
  suppressMessages(
    ggplot(data.frame(ate = ates), aes(x = ate)) +
      geom_histogram(bins = 30, colour = "black", fill = "skyblue") +
      geom_vline(xintercept = ci, linetype = "dashed") +
      labs(title = plot_title, x = "ATE", y = "Frequency") +
      theme_minimal(base_size = 14)
  )
  
  list(ates = ates, ci = ci, sd = sd_est)
}
