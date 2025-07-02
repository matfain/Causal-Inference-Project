# ============================================================================
#  CATE helpers – unified for T- and S-learners (RF / XGB)
#  (c) 2025  –  Matan · Hillel · Yael
# ============================================================================

library(tidyverse)
library(ranger)
library(xgboost)

# ---------------------------------------------------------------------------
# 0.  Helpers: align factor levels & build grids -----------------------------
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

build_grid <- function(df, modifier, n_points = 20) {
  x <- df[[modifier]]
  if (is.numeric(x))
    seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n_points)
  else
    levels(factor(x))
}

# ---------------------------------------------------------------------------
# 1.  Universal probabilistic learner wrapper (RF / XGB) --------------------
# ---------------------------------------------------------------------------

fit_prob_learner <- function(df, formula,
                             learner = c("rf", "xgb"),
                             # RF:
                             rf_ntree = 500,
                             # XGB:
                             xgb_nrounds = 300, xgb_eta = 0.05, xgb_depth = 6) {
  
  learner <- match.arg(learner)
  
  if (learner == "rf") {
    y_col <- all.vars(formula)[1]
    df_rf <- df
    if (is.numeric(df_rf[[y_col]]))
      df_rf[[y_col]] <- factor(df_rf[[y_col]], levels = c(0, 1))
    
    mod <- ranger::ranger(
      formula      = formula,
      data         = df_rf,
      probability  = TRUE,
      num.trees    = rf_ntree,
      mtry         = floor(sqrt(ncol(df_rf) - 2))
    )
    pred_prob <- function(newdata)
      predict(mod, newdata)$predictions[, "1"]
    
    return(list(fit = mod, pred_prob = pred_prob))
  }
  
  # ----------------------------  XGB  --------------------------------------
  X <- model.matrix(formula, data = df)[, -1, drop = FALSE]
  y <- df[[all.vars(formula)[1]]]; if (is.factor(y)) y <- as.numeric(y) - 1
  dtrain <- xgboost::xgb.DMatrix(data = X, label = y)
  
  mod <- xgboost::xgb.train(
    params = list(objective = "binary:logistic",
                  eta = xgb_eta, max_depth = xgb_depth,
                  eval_metric = "logloss"),
    data   = dtrain,
    nrounds = xgb_nrounds,
    verbose = 0)
  
  train_cols <- mod$feature_names        # cache once
  pred_prob  <- function(newdata) {
    newdata <- align_levels(newdata, df)
    M <- model.matrix(formula, data = newdata)[, -1, drop = FALSE]
    
    # pad / drop / reorder columns
    miss  <- setdiff(train_cols, colnames(M))
    if (length(miss))
      M <- cbind(M,
                 matrix(0, nrow(M), length(miss),
                        dimnames = list(NULL, miss)))
    extra <- setdiff(colnames(M), train_cols)
    if (length(extra))
      M <- M[, setdiff(colnames(M), extra), drop = FALSE]
    M <- M[, train_cols, drop = FALSE]
    
    stats::predict(mod, M)
  }
  
  list(fit = mod, pred_prob = pred_prob)
}

# ---------------------------------------------------------------------------
# 2.  T-learner --------------------------------------------------------------
# ---------------------------------------------------------------------------

fit_t_learner <- function(df, treatment, outcome, covars,
                          learner = c("rf", "xgb"), ...) {
  
  rhs  <- paste(c(treatment, covars), collapse = " + ")
  fmla <- as.formula(paste(outcome, "~", rhs))
  
  mod1 <- fit_prob_learner(df[df[[treatment]] == 1, ], fmla, learner, ...)
  mod0 <- fit_prob_learner(df[df[[treatment]] == 0, ], fmla, learner, ...)
  
  mu1  <- mod1$pred_prob(df)
  mu0  <- mod0$pred_prob(df)
  
  tibble(mu1 = mu1, mu0 = mu0, tau = mu1 - mu0,
         fit1 = list(mod1), fit0 = list(mod0))
}

# ---------------------------------------------------------------------------
# 3.  S-learner --------------------------------------------------------------
# ---------------------------------------------------------------------------

fit_s_learner <- function(df, treatment, outcome, covars,
                          learner = c("rf", "xgb"), ...) {
  
  rhs  <- paste(c(treatment, covars), collapse = " + ")
  fmla <- as.formula(paste(outcome, "~", rhs))
  
  base <- fit_prob_learner(df, fmla, learner, ...)
  
  df1 <- df; df1[[treatment]] <- 1
  df0 <- df; df0[[treatment]] <- 0
  
  mu1 <- base$pred_prob(df1)
  mu0 <- base$pred_prob(df0)
  
  tibble(mu1 = mu1, mu0 = mu0, tau = mu1 - mu0,
         base_model = list(base))
}

# ---------------------------------------------------------------------------
# 4.  Generic mean-CATE engine  ---------------------------------------------
# ---------------------------------------------------------------------------

cate_grid_mean_generic <- function(f_mu1, f_mu0, df, modifier, grid_vals) {
  map_dfr(grid_vals, function(g) {
    nd <- df; nd[[modifier]] <- g; nd <- align_levels(nd, df)
    tibble(mod_value = g, cate_mean = mean(f_mu1(nd) - f_mu0(nd)))
  })
}

cate_grid_mean <- function(fit0, fit1, df, treatment,
                           modifier, grid_vals) {
  cate_grid_mean_generic(fit1$pred_prob, fit0$pred_prob,
                         df, modifier, grid_vals)
}

cate_grid_mean_s <- function(base_model, df, treatment,
                             modifier, grid_vals) {
  cate_grid_mean_generic(
    function(nd) { r <- nd; r[[treatment]] <- 1; base_model$pred_prob(r) },
    function(nd) { r <- nd; r[[treatment]] <- 0; base_model$pred_prob(r) },
    df, modifier, grid_vals)
}

# ---------------------------------------------------------------------------
# 5.  Plot helpers (unchanged) ----------------------------------------------
# ---------------------------------------------------------------------------
plot_cate_vs_modifier <- function(df, cate_vec, modifier, learner_name, modifier_labels, ylim_range = c(-0.2, 0.2)) {
  
  plot_df <- df %>% mutate(tau = cate_vec, mod = .data[[modifier]])
  
  if (modifier == "sofa_first_cat") {
    plot_df$mod <- factor(plot_df$mod, levels = c("0-2.5", "2.5-5", "5-7.5", "7.5-10", "10+"))
  }
  
  x_label <- modifier_labels[[modifier]] %||% modifier
  
  if (is.numeric(plot_df$mod)) {
    p <- ggplot(plot_df, aes(mod, tau)) +
      geom_point(alpha = .2) +
      geom_smooth(method = "loess", se = FALSE, colour = "red") +
      ylim(ylim_range[1], ylim_range[2]) +
      labs(x = x_label, y = "Estimated CATE")
  } else {
    p <- ggplot(plot_df, aes(mod, tau, fill = mod)) +
      geom_boxplot(outlier.alpha = .15) +
      ylim(ylim_range[1], ylim_range[2]) +
      guides(fill = "none") +
      labs(x = x_label, y = "Estimated CATE")
  }
  print(p); invisible(p)
}

plot_cate_grid <- function(grid_df, modifier, learner_name) {
  
  if (is.numeric(grid_df$mod_value)) {
    p <- ggplot(grid_df, aes(mod_value, cate_mean)) +
      geom_line() + geom_point() +
      labs(title = sprintf("Mean CATE vs %s (%s)", modifier, learner_name),
           x = modifier, y = "Mean CATE")
  } else {
    p <- ggplot(grid_df, aes(mod_value, cate_mean)) +
      geom_col(fill = "steelblue") +
      labs(title = sprintf("Mean CATE by %s (%s)", modifier, learner_name),
           x = modifier, y = "Mean CATE")
  }
  print(p); invisible(p)
}
