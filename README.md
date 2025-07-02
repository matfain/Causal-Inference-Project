# 50 Shades of Zero – *the non-causal effect of IAC* 🎯

<!-- cool cover art -->
<div align="center">
  <img src="cover_image.png" width="450" alt="50 Shades of Zero">
</div>

&nbsp;

> **Causal question**  
> Does placing an indwelling arterial catheter (IAC) alter 28-day mortality in haemodynamically-stable, mechanically-ventilated ICU patients?

**Project Workflow**  
* 🧹 preprocessing of the **MIMIC-II** cohort dataset 
* 🎯 estimate the **ATE** with:  
  * ⚖️ IPTW using PS weighted difference in means & PS weighted outcome modeling  
  * 📏 Standardisation 
* 🧩 estimate **CATE** with metalearners:
  * 🌲 + 🦾 S-Learner using RF & XGBoost base learners
  * 🌲 + 🦾 T-Learner using RF & XGBoost base learners 


---

## 📂 Folder tree
```
Causal-Inference-Project
├── data
│   ├── full_cohort_data.csv
│   └── data_dictionary.txt
├── utils
│   ├── data_pipelines_utils.R
│   ├── plot_utils.R
│   ├── ate_utils_iptw.R
│   ├── ate_utils_standardization.R
│   └── cate_utils.R
├── ATE_IPTW_modeling.Rmd
├── ATE_standardization_modeling.Rmd
├── CATE_modeling_S_learner.Rmd
├── CATE_modeling_T_learner.Rmd
├── basic_eda.Rmd
├── requirements.R
└── README.md
```

---

## 🚀 How to run
1. **Clone**
   ```bash
   git clone https://github.com/matfain/Causal-Inference-Project.git
   cd Causal-Inference-Project
   ```

2. **Install packages**  
   ```r
   # in R / RStudio
   source("requirements.R")   # installs all dependencies required for the project
   ```

3. **Execute notebooks** – Open & Run the notebooks in the following order: 
   * `ATE_IPTW_modeling.Rmd` – IPTW pipeline for ATE estimation
   * `ATE_standardization_modeling.Rmd` – g-formula pipeline for ATE estimation  
   * `CATE_modeling_S_learner.Rmd` & `CATE_modeling_T_learner.Rmd` – Metalearners pipelines for CATE estimation

All notebooks auto-source helper functions from `utils/`, read the cohort from `data/`.
