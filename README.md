# Analysis of the Under-reported Antimicrobial Gene - *blaSCO-1*

This repository contains the data, scripts, and analysis results for the manuscript:

### **The Global Spread and Clinical Relevance of the Under-reported *blaSCO-1* Resistance Gene**



## Description

This project analyzes the temporal, host, and geographic distribution of the **blaSCO-1** β-lactamase gene using metadata from the NCBI Pathogen Detection Database. 

- The goal is to evaluate the clinical relevance and emergence patterns of *blaSCO-1* across bacterial species

> **Note:** Many isolates analyzed in this project have not been described in the literature. This analysis offers a first look into the global and cross-host spread of blaSCO-1 based on public metadata.



## Research Questions

- What is the global distribution of the *blaSCO-1* gene across bacterial species and hosts?
- How does its temporal distribution inform its clinical relevance?
- What host–species combinations are associated with earlier or more recent detection?



## Dataset Overview

- **Source:** [NCBI Pathogen Detection → MicroBIGG-E](https://www.ncbi.nlm.nih.gov/pathogens/microbigge/) using the keyword `"blasco-1"`
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



## Repository Structure

Here is the file tree structure of this repository:

```
Project
├── Blasco-1-Analysis.Rmd               # Top-level R Markdown file with integrated code and reporting
├── Code                                # Contains source scripts for analysis
│   └── blasco_analysis.R               # R script with codes only 
├── Data                                # Stores both raw and processed data
│   ├── blasco_cleaned_data.csv         # Processed and standardized data
│   └── blasco_rawdata.csv              # Original data
├── README.md                           # Top-level documentation with project overview and usage instructions
├── Results                             # All output files including plots, reports, and statistical summaries
│   ├── Blasco-1-Analysis_files         # Contains auto-generated gfm files from knitting
│   │   └── figure-gfm                  # Contains auto-generated gfm figures corresponding to chunk code
│   │       ├── unnamed-chunk-10-1.png  # Output from code chunk 10
│   │       ├── unnamed-chunk-3-1.png   # Output from code chunk 3
│   │       ├── unnamed-chunk-3-2.png   # -----
│   │       ├── unnamed-chunk-4-1.png
│   │       ├── unnamed-chunk-5-1.png
│   │       ├── unnamed-chunk-5-2.png
│   │       ├── unnamed-chunk-5-3.png
│   │       ├── unnamed-chunk-6-1.png
│   │       ├── unnamed-chunk-7-1.png
│   │       ├── unnamed-chunk-7-2.png
│   │       └── unnamed-chunk-8-1.png
│   ├── Documents                       # Final output reports in various formats
│   │   ├── Blasco-1-Analysis.html      # HTML report from knitting the Rmd file
│   │   └── Blasco-1-Analysis.md        # GFM markdown report for GitHub rendering
│   ├── Figures                         # Custom/exported figures manually saved (in JPEG format)
│   │   ├── HostSpecies_distribution_plots.jpg
│   │   ├── Host_isolate_source_association.jpg
│   │   ├── Host_species_interaction.jpg
│   │   ├── collection_date_significance.jpg
│   │   ├── global_distribution_map.jpg
│   │   ├── normality_plot.jpg
│   │   └── temporal_distribution_map.jpg
│   └── Statistics_Outputs               # CSV files containing test results and summary tables
│       ├── emmeans_species_host.csv
│       ├── pairwise_wilcox_Host.csv
│       └── summary_table.csv
└── blasco-1-analysis.Rproj             # Top-level RStudio project file to manage project settings

```

## How to Run the Analysis

**1. Clone the repository:**

   `git clone <https://github.com/NVI0001/blasco-1-analysis>`
   

**2. Install Dependencies**

The analysis was performed in R (version 4.4.1) and relies on several R packages. Make sure they are installed before running the code. You can install them using the command below:

```
install.packages(c(
  "stringr",      # String manipulation
  "dplyr",        # Data wrangling
  "tidyverse",    # Collection of core tidy data tools
  "ggplot2",      # Plotting system
  "forcats",      # Factor manipulation
  "maps",         # Map data (e.g., world outlines)
  "lme4",         # Mixed-effects models
  "ggpubr",       # Publication-ready visualizations
  "ggrepel",      # Non-overlapping text labels in ggplot
  "emmeans",      # Estimated marginal means
  "fs",           # File and directory handling
  "multcomp"      # Multiple comparisons and post-hoc tests
))

```


## Running the Analysis

To reproduce the results;

**Preferred Method:**

  Open Blasco-1-Analysis.Rmd in RStudio and click "Knit". This will create an HTML file that includes:

  - Data visualizations
  - Statistical results
  - Key findings

All analysis steps are compiled into one readable document.

**Alternative Method:**

  Open the blasco_analysis.R script (located in the Code/blasco_analysis.R) and run the code manually.
  > Note that: This script performs the analysis but does not include result interpretation as found in the .Rmd report.


## Results Summary

- The *blaSCO-1* antimicrobial resistance gene was detected in a wide range of bacterial species and geographic locations,
- The first detection was in the 1980s with a notable increase in reports after 2013. 
- Over 68% of detections were from human clinical isolates, highlighting its clinical significance. 
- *Klebsiella pneumoniae* was the first and most common bacterial species carrying the gene, followed by *E. coli*, *P. aeruginosa*, and *Salmonella enterica*. 
- The results indicate that *blaSCO-1* is more widespread than previously reported and call for expanded One Health surveillance.


## Citation

Please use the following DOI to cite this data and associated analyses/results.  

[![DOI](https://zenodo.org/badge/924392893.svg)](https://doi.org/10.5281/zenodo.14957149)

