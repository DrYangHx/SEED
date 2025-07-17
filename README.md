# SEED: A Novel Method for Identifying Cancer Prognosis Target Genes Without Control Samples

## 🔍 Project Overview

This repository provides the full analysis pipeline and datasets used in our study:  
**“SEED: A Novel Method for Identifying Cancer Prognosis Target Genes Without Control Samples.”**
The dataset can be downloaded from Baidu Netdisk:
Link: https://pan.baidu.com/s/1GuuHmOaitdvjyB3QHcd-Xw
Extraction Code: 1234
## 📁 Repository Structure

The GitHub repository contains **four core modules**:



# Clinical Variable Processing and Analysis


## 1. Clinical Variable Processing

The SEED method handles three types of clinical variables:

- **Ordered Categorical Variables**  
  Analyzed using Spearman's rank correlation, a non-parametric method that assesses monotonic relationships using rank statistics. This method is particularly suitable for ordered categorical variables (e.g., clinical stage) due to its robustness against outliers and non-normality, offering superior reliability compared to Pearson correlation for such data types.

- **Survival Data**  
  Analyzed using Cox regression, a semi-parametric model that accounts for survival time and event occurrence. This is the preferred method for analyzing time-to-event variables such as overall survival (OS) or progression-free survival (PFS), which also involve censoring considerations.

- **Unordered Categorical Variables**  
  Analyzed using Logistic regression, a probabilistic, nonlinear regression model designed specifically for categorical outcomes. It is the standard approach for analyzing unordered categorical dependent variables (e.g., objective response rate (ORR), viral positivity/negativity) by modeling binary or nominal outcomes through log-odds transformations.

## 2. Gene Expression Normalization

Scripts for preprocessing raw RNA-seq data, including normalization and missing value handling.

## 3. Feature Selection

Tools to perform:
- Statistical filtering (p-value threshold setting)
- Pathway enrichment analysis
- Candidate gene screening

## 4. Model Construction and Validation

Code for building prognostic models using:
- Lasso regression
- Multiple survival analysis methods (e.g., Cox regression, Kaplan-Meier, time-dependent ROC)


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
