#  A Multi-Study Modeling Approach

[![R-project](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![C++](https://img.shields.io/badge/Language-C++-orange.svg)](https://isocpp.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Official repository for the research paper **"Unveiling Mutational Signatures in Colon Adenocarcinoma: A Multi-Study Modeling Approach"** (2026).

This repository contains the Bayesian multi-study Non-Negative Matrix Factorization (NMF) MCMC sampler, dataset preprocessing code for TCGA-COAD/HNSC cohorts, and complete downstream clinical evaluation pipelines (Kaplan-Meier survival curves and Cox Forest plots).

---

## 🔬 Key Methodological Highlights

- **Multi-Study Framework**: Simultaneously models multiple genomic cohorts by decomposing mutation matrices into:
  - **Shared Signatures**: Universal biological processes across cohorts (e.g., aging/SBS1, POLE mutations/SBS10b).
  - **Study-Specific Signatures**: Cohort-exclusive features (e.g., SBS21 in homogeneous tumors; SBS103 and SBS46 in heterogeneous tumors).
- **Intratumoral Heterogeneity**: Stratification of TCGA-COAD tumors into **Homogeneous** and **Heterogeneous** subgroups based on the **MATH Score** (*Mutant-Allele Tumor Heterogeneity*).
- **High Performance**: C++ MCMC sampler integrated into R via `Rcpp` and `RcppArmadillo`.

---

## 📁 Repository Structure

```text
├── R/
│   ├── math_score.R         # Calculates MATH Score from TCGA MAF files
│   ├── cosmic_matching.R    # Cosine similarity calculation against COSMIC v3.6
│   ├── survival_analysis.R  # Kaplan-Meier survival analysis
│   └── forest_plots.R       # Cox proportional hazard ratios and Forest plots
├── src/
│   ├── mcmc_sampler.cpp     # High-performance C++ implementation of the MCMC sampler
│   └── RcppExports.cpp      # Auto-generated Rcpp bindings
├── data/
│   ├── readme.md            # Data access instructions (TCGA GDC Portal)
│   └── count_matrices.RData # Preprocessed mutation count matrices
├── figs/                    # Exported figures from paper (Boxplots, Heatmaps, Survival)
├── main.R                   # Master pipeline script execution
├── LICENSE                  # Open-source MIT License
└── README.md                # Project documentation
