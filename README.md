# Analysis of the Underreported Antimicrobial Gene - *blaSCO-1*

This repository contains the data, scripts, and analysis results for the manuscript:

### **The Global Spread and Clinical Relevance of the Underreported *blaSCO-1* Resistance Gene**



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



## File Structure
Here is the file tree structure of the repository:

```
.
├── Blasco-1-Analysis.Rmd               # Top-level R Markdown file with integrated code and reporting
├── Code                                # Contains source scripts for analysis
│   └── blasco_analysis.R               # R script with codes only 
├── Data                                # Stores both raw and processed data
│   ├── blasco_cleaned_data.csv
│   └── blasco_rawdata.csv
├── README.md                           # Top-level documentation with project overview and usage instructions
├── Results                             # All output files including plots, reports, and statistical summaries
│   ├── Blasco-1-Analysis_files     
│   │   └── figure-gfm                  # Contains auto-generated gfm figures
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

## Reproducing the Analysis

To run the full analysis, open `Blasco-1-Analysis.Rmd` in **RStudio** and knit it to **HTML**.

## Installation and Setup

1. Clone the repository:  
   `git clone <https://github.com/NVI0001/blasco-1-analysis>`
   
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

```


## Citation

Please use the following DOI to cite this data and associated functions/analyses.  
**DOI:** [Your DOI Here]
