# requirements.R

# List all packages used in the project
packages <- c(
  "tidyverse",
  "skimr",
  "tableone",
  "GGally",
  "naniar",
  "ggridges",
  "ggcorrplot",
  "factoextra",
  "vcd",
  "glmnet",
  "dplyr",
  "knitr",
  "tibble",
  "ranger",
  "xgboost",
  "doParallel",
  "foreach",
  "margins"
)

install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# install any missing package
invisible(lapply(packages, install_if_missing))
