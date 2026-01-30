# GR2D2: Grouped R2D2 shrinkage priors

This repository provides R code implementing grouped R2D2 prior, as proposed by the following paper.

Yanchenko, E., Irie, K. and Sugasawa, S. (2024). The Group R2D2 Shrinkage Prior for Sparse Linear Models with Grouped Covariates. [arXiv:2412.15293.](https://arxiv.org/abs/2412.15293) *Statistics and Computing* 36, article number 57.

The repository includes the following 7 files.

- `gR2D2.R`: Scriot to implement the proposed CSM model for graphical modeling
- `Competitors.R`: Scriot to implement existing shrinkage pirors
- `Example.R`: Scriot to run a oneshot simulation
- `Abalone-data.R`: Script for analysis of Abalone data
- `Abalone-data-validation.R`: Script for out-of-sample validation of Abalone data
- `Monte-Carlo-Simulation.R`: Script for Monte Carlo simulation in Yanchenko et al. (2024)
- `Summary-tale.R`: Script to summarize Monte Carlo simulation results 

