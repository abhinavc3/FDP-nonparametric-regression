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

- R (version used for this reproduction environment: `4.4.1`)
- R packages used in the repository:
  - `wavethresh` `4.7.3`
  - `ggplot2` `3.5.1`
  - `dplyr` `1.1.4`
  - `tidyr` `1.3.0`
  - `Metrics` `0.1.4`
  - `GGally` `2.2.1`
  - `RColorBrewer` `1.1-3`
  - `patchwork` `1.2.0`
  - `gridExtra` `2.3`
  - `foreign` `0.8-83`
  - `tidyverse` `2.0.0`
  - `haven` `2.5.4`
  - `xtable` `1.8-4`
  - `reshape2` `1.4.4`

---

## ▶️ Usage

To reproduce analyses:

1. Run the desired script directly from the checked-out repository; the analysis and simulation scripts are written to resolve helper files, inputs, and outputs relative to the script location rather than to a machine-specific working directory.
2. Use the scripts under `Data Analysis/` for data applications and the scripts under `Simulations/` for the main-paper and supplement experiments.

Example:

```bash
Rscript "Numerical Analysis for FED_REG/Simulations/main article sims/Figure 2 code/n_over_m_simulation.r"
