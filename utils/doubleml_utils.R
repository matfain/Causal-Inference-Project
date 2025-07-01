# ============================================================================
#  Double-ML / R-Learner helpers - FIXED VERSION
#  (K half-splits, default K = 5, 50-50 train/test each split)
# ============================================================================

library(tidyverse)     # pipes, dplyr, purrr
library(ranger)        # random-forest outcome / τ model
library(xgboost)       # XGB τ model
library(glmnet)        # (only for factor handling)

# ---------------------------------------------------------------------------
# 1  Make K random 50-50 splits  --------------------------------------------
# ---------------------------------------------------------------------------
make_half_splits <- function(n, K = 5, seed = 1) {
  set.seed(seed)
  map(seq_len(K), ~ sample(seq_len(n), size = floor(n / 2)))
}

# ---------------------------------------------------------------------------
# 2  Outcome model  m(x)  (RF) ----------------------------------------------
# ---------------------------------------------------------------------------
fit_outcome_rf <- function(df, fmla, ntree = 500) {
  y_name <- all.vars(fmla)[1]
  
  df2 <- df
  if (is.numeric(df2[[y_name]]))
    df2[[y_name]] <- factor(df2[[y_name]], levels = c(0, 1))
  
  mod <- ranger::ranger(
    formula     = fmla,
    data        = df2,
    probability = TRUE,
    num.trees   = ntree,
    mtry        = floor(sqrt(ncol(df2) - 2))
  )
  list(
    fit  = mod,
    pred = function(newdata) predict(mod, newdata)$predictions[, "1"]
  )
}

# ---------------------------------------------------------------------------
# 3  Propensity model  p(x)  (plain logistic glm) ---------------------------
# ---------------------------------------------------------------------------
fit_propensity_glm <- function(df, fmla) {
  mod <- glm(fmla, data = df, family = binomial)
  list(
    fit  = mod,
    pred = function(newdata) stats::predict(mod, newdata, type = "response")
  )
}

# ---------------------------------------------------------------------------
# 4  Pseudo outcome & weights - FIXED VERSION -------------------------------
#     ψ = (Y - m̂) / (A - p̂) ;   w = (A - p̂)^2
# ---------------------------------------------------------------------------
make_pseudo <- function(df_hold, m_hat, p_hat, outcome, treatment, 
                        min_propensity = 0.01, max_propensity = 0.99) {
  
  y  <- df_hold[[outcome]]
  a  <- df_hold[[treatment]]
  m  <- m_hat(df_hold)
  p  <- p_hat(df_hold)
  
  # Clip propensity scores to avoid division by zero
  p <- pmax(min_propensity, pmin(max_propensity, p))
  
  # Check for missing values in predictions
  if (any(is.na(m)) || any(is.na(p))) {
    warning("Missing values detected in nuisance predictions")
  }
  
  # Calculate pseudo outcomes
  denominator <- a - p
  psi <- (y - m) / denominator
  w   <- denominator^2
  
  
  # Remove observations with invalid pseudo outcomes
  valid_idx <- is.finite(psi) & is.finite(w) & !is.na(psi) & !is.na(w)
  
  if (sum(valid_idx) < length(psi)) {
    warning(paste("Removed", sum(!valid_idx), "observations with invalid pseudo outcomes"))
  }
  
  tibble(
    psi = ifelse(valid_idx, psi, NA),
    w = ifelse(valid_idx, w, NA),
    valid = valid_idx
  )
}

# ---------------------------------------------------------------------------
# 5  Align factor/character levels  -----------------------------------------
# ---------------------------------------------------------------------------
align_levels <- function(newdata, reference_df) {
  for (nm in names(newdata)) {
    if (is.factor(reference_df[[nm]])) {
      newdata[[nm]] <- factor(newdata[[nm]],
                              levels = levels(reference_df[[nm]]))
    } else if (is.character(reference_df[[nm]])) {
      newdata[[nm]] <- factor(newdata[[nm]],
                              levels = unique(reference_df[[nm]]))
    }
  }
  newdata
}

# ---------------------------------------------------------------------------
# 6  τ-model learners  (RF  /  XGB) - FIXED VERSION ------------------------
# ---------------------------------------------------------------------------
fit_tau_model <- function(df, psi, w, covars,
                          learner = c("rf", "xgb"),
                          rf_ntree = 500,
                          xgb_nrounds = 300, xgb_eta = 0.05, xgb_depth = 6) {
  
  learner <- match.arg(learner)
  
  # Remove observations with missing pseudo outcomes
  valid_idx <- !is.na(psi) & !is.na(w) & is.finite(psi) & is.finite(w)
  
  if (sum(valid_idx) == 0) {
    stop("No valid observations for tau model fitting")
  }
  
  if (sum(valid_idx) < length(psi)) {
    # print the invalid values and their indices
    invalid_indices <- which(!valid_idx)
    invalid_values <- psi[!valid_idx]
    warning(paste("Invalid pseudo outcomes at indices:", 
                  paste(invalid_indices, collapse = ", "),
                  "with values:", 
                  paste(invalid_values, collapse = ", ")))
    warning(paste("Using", sum(valid_idx), "out of", length(psi), "observations for tau model"))
  }
  
  df_aug  <- df[valid_idx, ] %>% 
    mutate(psi = psi[valid_idx], weight = w[valid_idx])
  
  if (learner == "rf") {
    fmla <- as.formula(paste("psi ~", paste(covars, collapse = " + ")))
    mod  <- ranger::ranger(
      formula     = fmla,
      data        = df_aug,
      case.weights= df_aug$weight,
      num.trees   = rf_ntree
    )
    return(list(
      fit  = mod,
      pred = function(newdata) predict(mod, newdata)$predictions
    ))
  }
  
  # ------- XGB -------------------------------------------------------------
  X <- model.matrix(
    ~ . -1,
    data = align_levels(df_aug[, covars, drop = FALSE], df_aug)
  )
  dtrain <- xgboost::xgb.DMatrix(data = X, label = df_aug$psi,
                                 weight = df_aug$weight)
  
  mod <- xgboost::xgb.train(
    params = list(
      objective   = "reg:squarederror",
      eta         = xgb_eta,
      max_depth   = xgb_depth,
      subsample   = 0.7,
      colsample_bytree = 0.7
    ),
    data     = dtrain,
    nrounds  = xgb_nrounds,
    verbose  = 0
  )
  list(
    fit  = mod,
    feature_names = colnames(X),              # keep for prediction
    pred = function(newdata) {
      newX <- model.matrix(~ . -1,
                           data = align_levels(newdata[, covars, drop = FALSE], df_aug))
      # column reconcile
      miss <- setdiff(mod$feature_names, colnames(newX))
      if (length(miss) > 0)
        newX <- cbind(newX,
                      matrix(0, nrow(newX), length(miss),
                             dimnames = list(NULL, miss)))
      newX <- newX[, mod$feature_names, drop = FALSE]
      stats::predict(mod, newX)
    }
  )
}

# ---------------------------------------------------------------------------
# 7  Full R-learner pipeline - FIXED VERSION -------------------------------
# ---------------------------------------------------------------------------
r_learner <- function(df, treatment, outcome, covars,
                      K = 5,
                      learners = list(outcome    = "rf",
                                      propensity = "glm",
                                      tau        = "rf"),
                      seed = 1,
                      min_propensity = 0.01,
                      max_propensity = 0.99) {
  
  n          <- nrow(df)
  fold_index <- make_half_splits(n, K, seed)
  
  fmla_y  <- as.formula(paste(outcome,   "~", paste(covars, collapse = " + ")))
  fmla_a  <- as.formula(paste(treatment, "~", paste(covars, collapse = " + ")))
  
  # store pseudo-outcomes from each split
  psi_all <- rep(NA_real_, n)
  w_all   <- rep(NA_real_, n)
  
  for (k in seq_len(K)) {
    
    test_id  <- fold_index[[k]]
    train_id <- setdiff(seq_len(n), test_id)
    
    df_train <- df[train_id, , drop = FALSE]
    df_test  <- df[test_id,  , drop = FALSE]
    
    # ---- nuisances --------------------------------------------------------
    m_hat <- fit_outcome_rf(df_train, fmla_y)$pred
    p_hat <- fit_propensity_glm(df_train, fmla_a)$pred
    
    # ---- pseudo outcomes on hold-out -------------------------------------
    pseudo <- make_pseudo(df_test, m_hat, p_hat,
                          outcome   = outcome,
                          treatment = treatment,
                          min_propensity = min_propensity,
                          max_propensity = max_propensity)
    
    psi_all[test_id] <- pseudo$psi
    w_all[test_id]   <- pseudo$w
  }
  
  # Check how many valid pseudo outcomes we have
  valid_psi <- sum(!is.na(psi_all))
  cat("Valid pseudo outcomes:", valid_psi, "out of", n, "\n")
  
  if (valid_psi < n * 0.5) {
    warning("Less than 50% of observations have valid pseudo outcomes")
  }
  
  # ---- τ-model on pseudo data with valid observations ------------------
  tau_fit <- fit_tau_model(df, psi_all, w_all, covars,
                           learner = learners$tau)
  
  tau_vec <- tau_fit$pred(df)
  
  list(tau_vec   = tau_vec,
       tau_model = tau_fit,
       psi       = psi_all,
       w         = w_all,
       folds     = fold_index,
       valid_obs = sum(!is.na(psi_all)))
}

# ---------------------------------------------------------------------------
# 8  Utility: mean-CATE over grid using τ-model ------------------------------
# ---------------------------------------------------------------------------
cate_grid_mean_tau <- function(tau_fit, df, modifier, grid_vals) {
  
  map_dfr(grid_vals, function(g) {
    
    newdata <- df
    newdata[[modifier]] <- g
    newdata <- align_levels(newdata, df)
    
    tibble(mod_value = g,
           cate_mean = mean(tau_fit$pred(newdata)))
  })
}