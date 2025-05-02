# requirements.R

packages <- c(
  "tidyverse",
  "skimr",
  "tableone",
  "GGally",
  "naniar",
  "ggridges",
  "ggcorrplot",
  "factoextra",
  "vcd"
)

install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

invisible(lapply(packages, install_if_missing))
