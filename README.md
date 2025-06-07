# Causal-Inference-Project

This repository contains code and analyses for a university project on **Causal Inference**, primarily using **R** and **R Markdown**.  
The goal is to explore, estimate, and interpret causal effects using modern statistical and machine learning methods.

## 📊 Dataset Overview

### Data Source
This project uses data from the **Multiparameter Intelligent Monitoring in Intensive Care II (MIMIC-II)** database, containing detailed clinical data from over 24,000 patients admitted to Beth Israel Deaconess Medical Center ICUs between 2001 and 2008.

### Study Population
From the full MIMIC-II database, we analyze **1,776 patients** who meet specific criteria:
- Adult patients requiring mechanical ventilation within 24 hours of ICU admission
- Mechanical ventilation lasting at least 24 hours  
- Hemodynamically stable (no vasopressor support needed)
- No sepsis diagnosis

### Treatment Variable
**Indwelling Arterial Catheter (IAC) placement** - defined as insertion of an invasive arterial catheter after initiation of mechanical ventilation. In our cohort, 44.6% (792 patients) received an IAC.

### Key Variables
The dataset includes:
- Patient demographics (age, gender, BMI)
- Clinical scores (SOFA, SAPS I)
- Comorbidities (CHF, COPD, stroke, etc.)
- Laboratory values (blood gases, chemistry panel)
- Outcomes (28-day mortality, ICU/hospital length of stay)
- Treatment indicator (aline_flg)

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
