# ============================================================================
#  CATE helpers – T-learner
# ============================================================================

# ---------------------------------------------------------------------------
# 1.  Fit a single base learner  --------------------------------------------
# ---------------------------------------------------------------------------
fit_base_model <- function(df, formula,
                           learner = c("rf", "xgb"),
                           rf_ntree = 500,
                           xgb_nrounds = 200, xgb_eta = 0.1, xgb_depth = 6) {
  
  learner <- match.arg(learner)
  
  if (learner == "rf") {
    outcome_name <- all.vars(formula)[1]
    
    df_rf <- df
    if (is.numeric(df_rf[[outcome_name]])) {
      df_rf[[outcome_name]] <- factor(df_rf[[outcome_name]], levels = c(0, 1))
    }
    
    mod <- ranger::ranger(
      formula      = formula,
      data         = df_rf,
      probability  = TRUE,
      num.trees    = rf_ntree,
      mtry         = floor(sqrt(ncol(df_rf) - 2))
    )
    
    predict_prob <- function(newdata) {
      pred <- predict(mod, newdata)$predictions
      # pred is an N×2 matrix; take column "1"
      pred[, "1"]
    }
    
    return(list(fit = mod, predict_prob = predict_prob))
  }
  
  
  if (learner == "xgb") {
    
    # --- design matrix (drop intercept) ----
    X <- model.matrix(formula, data = df)[, -1, drop = FALSE]
    y <- df[[all.vars(formula)[1]]]
    if (is.factor(y)) y <- as.numeric(y) - 1
    
    dtrain <- xgboost::xgb.DMatrix(data = X, label = y)
    
    mod <- xgboost::xgb.train(
      params   = list(
        objective   = "binary:logistic",
        eta         = xgb_eta,
        max_depth   = xgb_depth,
        eval_metric = "logloss"
      ),
      data     = dtrain,
      nrounds  = xgb_nrounds,
      verbose  = 0
    )
    
    # predict_prob <- function(newdata) {
    #   newX <- model.matrix(formula, data = newdata)[, -1, drop = FALSE]
    #   stats::predict(mod, newX)          # vector of P(Y=1)
    # }
    
    predict_prob <- function(newdata) {
      
      # 1. align factor levels
      newdata <- align_levels(newdata, df)
      
      # 2. build design matrix
      newX <- model.matrix(formula, data = newdata)[, -1, drop = FALSE]
      
      # 3. reconcile columns with training matrix
      train_cols <- mod$feature_names
      
      if (!identical(colnames(newX), train_cols)) {
        # add missing
        missing <- setdiff(train_cols, colnames(newX))
        if (length(missing) > 0) {
          newX <- cbind(newX,
                        matrix(0, nrow = nrow(newX), ncol = length(missing),
                               dimnames = list(NULL, missing)))
        }
        # drop extra
        extra <- setdiff(colnames(newX), train_cols)
        if (length(extra) > 0) {
          newX <- newX[, setdiff(colnames(newX), extra), drop = FALSE]
        }
        # reorder
        newX <- newX[, train_cols, drop = FALSE]
      }
      
      stats::predict(mod, newX)
    }
    
    return(list(fit = mod, predict_prob = predict_prob))
  }
  
  
  stop("Unsupported learner")
}

# ---------------------------------------------------------------------------
# 2.  Fit T-learner & return CATE vector  ------------------------------------
# ---------------------------------------------------------------------------
fit_t_learner <- function(df, treatment, outcome, covars, learner = "rf", ...) {
  
  # formula uses outcome as response; covars + all else in df
  rhs  <- paste(c(treatment, covars), collapse = " + ")
  fmla <- as.formula(paste(outcome, "~", rhs))
  
  # split
  df_1 <- df[df[[treatment]] == 1, ]
  df_0 <- df[df[[treatment]] == 0, ]
  
  # train models
  mod1 <- fit_base_model(df_1, fmla, learner, ... )
  mod0 <- fit_base_model(df_0, fmla, learner, ... )
  
  # predict mu1, mu0 **for every row**
  mu1 <- mod1$predict_prob(df)
  mu0 <- mod0$predict_prob(df)
  
  tibble(mu1 = mu1, mu0 = mu0, tau = mu1 - mu0,
         fit1 = list(mod1), fit0 = list(mod0))
}

# ---------------------------------------------------------------------------
# 3.  Plot CATE vs. effect modifier  -----------------------------------------
# ---------------------------------------------------------------------------
plot_cate_vs_modifier <- function(df, cate_vec, modifier, learner_name,
                                  bins = 10) {
  
  stopifnot(modifier %in% names(df))
  plot_df <- df %>% mutate(tau = cate_vec, mod = .data[[modifier]])
  
  if (is.numeric(plot_df$mod)) {
    # numeric modifier: bin + loess
    plot_df <- plot_df %>%
      mutate(bin = cut(mod, breaks = bins, include.lowest = TRUE))
    
    p <- ggplot(plot_df, aes(x = mod, y = tau)) +
      geom_point(alpha = 0.2) +
      geom_smooth(method = "loess", se = FALSE, colour = "red") +
      labs(
        title = sprintf("CATE vs %s  (%s)", modifier, learner_name),
        x = modifier, y = "Estimated CATE"
      )
    
  } else {
    # categorical modifier: boxplot
    p <- ggplot(plot_df, aes(x = mod, y = tau, fill = mod)) +
      geom_boxplot(outlier.alpha = 0.15) +
      guides(fill = "none") +
      labs(
        title = sprintf("CATE by %s  (%s)", modifier, learner_name),
        x = modifier, y = "Estimated CATE"
      )
  }
  
  print(p)
  invisible(p)
}

# ---------------------------------------------------------------------------
# 4.  Build grid for a modifier
# ---------------------------------------------------------------------------
build_grid <- function(df, modifier, n_points = 20) {
  x <- df[[modifier]]
  if (is.numeric(x)) {
    seq(min(x, na.rm = TRUE),
        max(x, na.rm = TRUE),
        length.out = n_points)
  } else {
    levels(factor(x))
  }
}

# ---- ensure new data carries the full factor/character level set -----------
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
# 5.  Mean-CATE over all patients for each grid value
# ---------------------------------------------------------------------------
cate_grid_mean <- function(fit0, fit1, df, treatment,
                           modifier, grid_vals) {
  
  purrr::map_dfr(grid_vals, function(g) {
    newdata <- df
    newdata[[modifier]] <- g
    
    newdata <- align_levels(newdata, df)
    
    # keep factor levels consistent
    if (is.factor(df[[modifier]])) {
      newdata[[modifier]] <- factor(newdata[[modifier]],
                                    levels = levels(df[[modifier]]))
    }
    
    mu0 <- fit0$predict_prob(newdata)
    mu1 <- fit1$predict_prob(newdata)
    
    tibble(mod_value = g, cate_mean = mean(mu1 - mu0))
  })
}

# ---------------------------------------------------------------------------
# 6.  Plot CATE grid
# ---------------------------------------------------------------------------
plot_cate_grid <- function(grid_df, modifier, learner_name) {
  
  if (is.numeric(grid_df$mod_value)) {
    p <- ggplot(grid_df, aes(x = mod_value, y = cate_mean)) +
      geom_line() + geom_point() +
      labs(title = sprintf("Mean CATE vs %s  (%s)", modifier, learner_name),
           x = modifier, y = "Mean CATE")
  } else {
    p <- ggplot(grid_df, aes(x = mod_value, y = cate_mean)) +
      geom_col(fill = "steelblue") +
      labs(title = sprintf("Mean CATE by %s  (%s)", modifier, learner_name),
           x = modifier, y = "Mean CATE")
  }
  print(p)
  invisible(p)
}
