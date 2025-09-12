# FDP Nonparametric Regression

This repository contains code and data used for numerical analysis and simulations in the context of federated nonparametric regression under privacy constraints.

---

## 📂 Repository Structure

- **Numerical Analysis for FED_REG/**
  - **Data Analysis/**
    - **heart dataset/**  
      Contains R scripts for loading, preprocessing, exploratory analysis, and estimation on the heart disease dataset. Includes raw data files in the `heart+disease/` subfolder with costs and metadata.
      - `wavelet_helper_functions.r` — helper functions for wavelet-based estimators  
      - `federated_est_heart.R` — federated estimation on heart data  
      - `EDA-heart.R` — exploratory data analysis  
      - `individual_est_heart.R` — individual (non-federated) estimation  
      - `load_data.R` — dataset loading script  
      - `heart.R` — main script for running heart data analysis  

    - **expiratory dataset/**  
      Scripts and datasets for lung function (FEV1) analysis. Includes R scripts and multiple PDF plots.  
      - `analysis.R` — analysis pipeline  
      - Several `.XPT` files — raw dataset inputs  
      - PDF files — plots for smoking vs non-smoking comparisons and wavelet estimations  

    - **homocysteine dataset/**  
      Analysis scripts and plots for homocysteine data with covariates like folate and B12.  
      - `analysis.R` — main analysis script  
      - `.XPT` files — raw dataset inputs  
      - Multiple PDFs — estimation results and plots  

  - **Simulations/**
    - **main article sims/**  
      Code and figures corresponding to the main article experiments.  
      - `Figure 1 code/` — scripts for generating Figure 1  
      - `Figure 2 code/` — scripts for generating Figure 2 (e.g., `n_over_m_simulation.r`, `eps_sim_logscale.r`)  
      - `Figures/` — output figures (PDFs)  

    - **additional sims for supplement/**  
      Supplementary simulations for robustness and tables.  
      - `table generation/` — scripts to generate LaTeX tables of simulation results  
      - `robustness sims/` — robustness analysis scripts (e.g., noise level, smoothness)  
      - `Figures/` — robustness result plots  

- **.git/**  
  Version control history and metadata.

---

## ⚙️ Requirements

- R (≥ 4.0 recommended)  
- R packages: `wavelets`, `ggplot2`, `data.table`, `foreign`, and others depending on specific scripts.

*(You may need to install additional packages as indicated in individual scripts.)*

---

## ▶️ Usage

To reproduce analyses:

1. Navigate to the dataset of interest inside `Data Analysis/` and run the corresponding `.R` scripts.
2. For simulation studies, go to `Simulations/` and run the scripts under `main article sims/` or `additional sims for supplement/`.

Example:

```bash
Rscript "Numerical Analysis for FED_REG/Simulations/main article sims/Figure 2 code/n_over_m_simulation.r"
