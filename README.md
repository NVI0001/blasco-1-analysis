# Analysis of the Underreported Antimicrobial Gene - blaSCO-1

This repository contains the data, scripts, and analysis results for the manuscript:

### **The Global Spread and Clinical Relevance of the Underreported blaSCO-1 Resistance Gene**

---

## Description

This project analyzes the temporal, host, and geographic distribution of the **blaSCO-1** β-lactamase gene using metadata from the NCBI Pathogen Detection Database. 

- The goal is to evaluate the clinical relevance and emergence patterns of blaSCO-1 across bacterial species

> **Note:** Many isolates analyzed in this project have not been described in the literature. This analysis offers a first look into the global and cross-host spread of blaSCO-1 based on public metadata.

---

## Research Questions

- What is the global distribution of the *blaSCO-1* gene across bacterial species and hosts?
- How does its temporal distribution inform its clinical relevance?
- What host–species combinations are associated with earlier or more recent detection?

---

## Dataset Overview

- **Source:** NCBI Pathogen Detection → MicroBIGG-E
- **Date Accessed:** 2024-04-08

### Metadata Fields:
- Scientific name  
- Host  
- Geographic origin  
- Isolation source & type  
- Collection date  
- Strand orientation  
- Method of detection  
- Element name: *blaSCO-1*

This metadata was downloaded in CSV format.

---
  


## File Tree
Here is the file tree structure of the repository:


## Reproducing the Analysis

To run the full analysis, open `Blasco-1-Analysis.Rmd` in **RStudio** and knit it to **HTML**.

## Installation and Setup

1. Clone the repository:  
   `git clone <repository-url>`
   
2. Install necessary R packages (if not already installed):
To run the analysis, the following R packages are required:
- `ggplot2`: For data visualization.
- `dplyr`: For data manipulation.
- `tidyr`: For data tidying.
- `vegan`: For ecological analyses (e.g., distance matrices).
- `knitr`: For generating reports from RMarkdown.
- `readr`: For reading CSV files.

You can install these packages by running the following in R:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "vegan", "knitr", "readr"))




## Citation
Please use the following DOI to cite this data and associated functions/analyses.  
**DOI:** [Your DOI Here]
