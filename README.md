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
├── R/
│   ├── Funcoes_analises_MEstudo.r       # Funções auxiliares e análises do estudo
│   └── De_novo.R                        # Script para extração de assinaturas de novo
├── src/
│   └── Gibbs_Sampler_D.cpp              # Implementação do amostrador de Gibbs em C++
├── data/
│   ├── all.tcga-clinical-indexed.tsv    # Dados clínicos indexados do TCGA
│   ├── COSMIC_SBS96_hg38_ordered.txt    # Referência mutacional COSMIC (SBS96 hg38)
│   ├── countsMATH_Homogeneos_COAD.txt   # Matrizes de contagem - COAD Homogêneos
│   └── countsMATH_Heterogeneos_COAD.txt # Matrizes de contagem - COAD Heterogêneos
├── figs/                                # Pasta vazia para receber futuros gráficos/outputs
├── Input_Multi_Estudo.R                 # Script principal (Atua como o "main.R" do projeto)
└── README.md                            # Documentação principal
