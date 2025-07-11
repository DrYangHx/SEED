# SEED: A Novel Method for Identifying Cancer Prognosis Target Genes Without Control Samples

## 🔍 Project Overview

This repository provides the full analysis pipeline and datasets used in our study:  
**“SEED: A Novel Method for Identifying Cancer Prognosis Target Genes Without Control Samples.”**

## 📁 Repository Structure

The GitHub repository contains **three core modules**:

1. **Gene Expression Normalization**  
   Scripts for preprocessing raw RNA-seq data, including normalization and missing value handling.

2. **Feature Selection**  
   Tools to perform statistical filtering (p-value threshold setting), pathway enrichment analysis, and candidate gene screening.

3. **Model Construction and Validation**  
   Code for building prognostic models using Lasso regression and multiple survival analysis methods (e.g., Cox, Kaplan-Meier, time-dependent ROC).

> 📌 For complete reproducibility, intermediate data and scripts are also available at:  
> - GitHub: https://github.com/DrYangHx/SEED  
> - Baidu Cloud: [insert link]  
> (Example: Full feature selection trace for the UCEC dataset from raw input to the final 8-gene model.)

## 📊 Parameter Documentation

All dataset-specific parameters used during the model-building process are provided in `result.xlsx`, including:

- p-value threshold for preliminary screening  
- q-value threshold in pathway enrichment  
- random seed for Lasso regression  
- the number of features retained after each step for each dataset  

## 💻 Software Environment

The SEED web application was developed in **R 4.3.3**, using the following major packages:

- `survival` (v3.7.0): for survival analysis  
- `glmnet` (v4.1.8): for Lasso regression  
- `ggplot2` (v3.5.1): for data visualization
