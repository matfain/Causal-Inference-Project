# Causal-Inference-Project

This repository contains code and analyses for a university project on **Causal Inference**, primarily using **R** and **R Markdown**.  
The goal is to explore, estimate, and interpret causal effects using modern statistical and machine learning methods.

## 📁 Project Structure

```
Causal-Inference-Project/
├── data/ # Raw data files
│ ├── data_dictionary.txt # Description of dataset variables
│ └── full_cohort_data.csv # Main cohort dataset for analysis
├── sandbox/ # Experimental notebooks and scratch work
│ └── base_eda_matan.Rmd # Initial exploratory data analysis
├── requirements.R # R packages required for the project
├── .gitignore # Ignore rules for Git version control
├── README.md # Project overview and structure
```

## 🛠 Technologies

- R (version 4.2 or above)  
- RStudio  
- R Markdown  
- Git & GitHub  

## 📦 Required R Packages

To install all required packages, run in your R console:

```r
source("requirements.R")
```

This script installs:
- tidyverse
- skimr
- tableone
- GGally
- naniar
- ggridges
- ggcorrplot
- factoextra
- vcd
- glmnet
- dplyr
- knitr
- tibble

## 📌 Notes
- Keep large or sensitive data outside Git (`data/` is tracked only for small public files)
- Use the `sandbox/` folder for temporary or in-progress analyses
