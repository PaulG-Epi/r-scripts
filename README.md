# Please read for Microbiome and Epigenetic Clock projects 

## **Microbiome Project:**
### **Nasal Cavity Microbial Makeup and the Influence on Psychiatric Symptoms Following Fire Exposure in Firefighters**


This repository contains R scripts for a study examining associations between the nasal microbiome and psychiatric outcomes (e.g., depression, anxiety, PTSD) in firefighters. The project evaluates how microbial composition relates to mental health measures using observational data.
Analyses include logistic regression, linear regression, and linear mixed-effects models to assess relationships between microbial features and psychiatric outcomes.

#### How to use:

1) install the following packages: openxlsx, writexl, plyr, aod, ggplot2, readxl, gtsummary, tidyr, lme4, lmerTest, gt, patchwork, flextable, officer, multtest, phyloseq, vegan, viridis
2) Download the script: Firefighter mental health microbiome code
3) Run the script in R studio

#### Code Overview

The code contains:
- data cleaning/organizing
- descriptive analyses 
- various regression modeling steps (i.e., logistic, linear, and linear mixed effects model steps)
- outputting tables and figures




## **Epigenetic Clock Project:**
### **Associations between occupational polycyclic aromatic hydrocarbon (PAH) exposure and epigenetic aging biomarkers in firefighters**


This repository contains R scripts for a study examining associations between environmental exposures (e.g., PAHs) and epigenetic age acceleration in firefighters using multiple DNA methylation clocks. The project evaluates how exposure measures relate to biological aging metrics.
Analyses include data cleaning, linear regression models, and leave-one-out approaches to assess the robustness of exposure–outcome relationships.

#### How to use:

1) install the following packages: readxl, dplyr, gtsummary, officer, flextable, tidyr, stringr, ggplot2, patchwork, tidyverse, readr, broom, purrr, DunedinPACE, methylclock, sesame, lme4, lmerTest, forcats
2) Download the script: MASTER_EPICLOCK_CODE
3) Run the script in R studio

#### Code Overview

The code contains:
- data cleaning/organizing
- descriptive analyses/tests
- epigenetic clock age estimates
- residual epigenetic age estimates
- linear regression modelling steps
- leave-one-out analyses
- outputting steps for tables and figures
