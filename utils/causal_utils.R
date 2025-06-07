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