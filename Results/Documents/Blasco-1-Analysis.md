## 1. Introduction

The *blaSCO-1* gene, first reported in *Salmonella* in the USA (Iduu et
al. 2024), has limited global documentation, with most reports
concentrated in a few regions. Its role in conferring resistance to
β-lactam antibiotics, including penicillins and cephalosporins, raises
concerns about its potential clinical impact. Despite its identification
in a few Enterobacteriaceae species like *Klebsiella* and *E. coli*, the
global prevalence and genetic context of *blaSCO-1* remain largely
underexplored.

- The goal is to evaluate the clinical relevance and emergence patterns
  of blaSCO-1 across bacterial species

### 1.1 Research Questions

- What is the global distribution of the *blaSCO-1* gene across
  bacterial species and hosts?
- How does its temporal distribution inform its clinical relevance?
- What host–species combinations are associated with earlier or more
  recent detection?

## 2. Materials and Methods

### 2.1 Data Collection

Metadata on bacterial isolates harboring the *blaSCO-1* gene were
retrieved from the NCBI Pathogen Detection Database, specifically
through the Pathogen Detection Microbial Browser for Identification of
Genetic and Genomic Elements (MicroBIGG-E) (accessed on 2024-04-08). The
dataset contains isolate-level information including scientific name,
host, geographic origin, isolation source, isolation type, collection
date, strand orientation, and method of detection. The metadata were
downloaded in CSV format and subsequently used for analysis.

### 2.2 Data Import and PreProcessing

Data were imported into R studio (v4.4..3) using the read.csv()
function, with missing values properly interpreted using na.strings =
c(““,”na”) and stringsAsFactors = TRUE to automatically convert
character columns to factors. This facilitates categorical data analysis
and visualization. Missing values were interpreted using na.strings =
c(““,”na”), and data were cleaned to ensure standardization.

``` r
# Load required libraries
library(stringr)      # For string manipulation functions (e.g., pattern matching, replacement, extraction)
library(dplyr)        # For data manipulation (e.g., filter, mutate, select, summarise, etc.)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyverse)    # Loads a collection of tidy data tools (includes dplyr, ggplot2, tidyr, readr, etc.)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ readr     2.1.5
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggplot2)      # For creating data visualizations and plots
library(forcats)      # For handling and reordering categorical/factor variables
library(maps)         # For drawing geographical maps (includes country/state boundaries)
```

    ## 
    ## Attaching package: 'maps'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     map

``` r
library(lme4)         # For fitting linear and generalized linear mixed-effects models
```

    ## Loading required package: Matrix
    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
library(ggpubr)       # For publication-ready plots, adds enhancements to ggplot2
library(ggrepel)      # For better label placement in ggplot2 (avoids label overlap)
library(emmeans)      # For estimating marginal means (EMMs) aka least-squares means from model outputs
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

``` r
library(fs)           # For generating file tree
library(multcomp)     # For multiple comparison procedures and simultaneous inference
```

    ## Loading required package: mvtnorm
    ## Loading required package: survival
    ## Loading required package: TH.data
    ## Loading required package: MASS
    ## 
    ## Attaching package: 'MASS'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select
    ## 
    ## 
    ## Attaching package: 'TH.data'
    ## 
    ## The following object is masked from 'package:MASS':
    ## 
    ##     geyser

``` r
#load data 
data <- read.csv("Data/blasco_rawdata.csv", na.strings = c("", "na"), stringsAsFactors = TRUE)

str(data) # View data structure
```

    ## 'data.frame':    768 obs. of  11 variables:
    ##  $ Scientific_name : Factor w/ 21 levels "Acinetobacter baumannii",..: 7 19 19 4 4 19 17 17 17 17 ...
    ##  $ Isolation._type : Factor w/ 3 levels "agricultural",..: 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ Location        : Factor w/ 46 levels "Argentina","Australia",..: 1 13 13 NA NA 13 46 46 46 46 ...
    ##  $ Isolation_source: Factor w/ 42 levels " wound","bile",..: 8 19 19 NA NA NA NA NA NA NA ...
    ##  $ Host            : Factor w/ 5 levels "chicken","dog",..: 4 4 4 NA NA 4 3 3 3 3 ...
    ##  $ Collection_date : int  1980 1988 1988 1988 1988 1988 1989 1989 1991 1991 ...
    ##  $ Isolate         : Factor w/ 506 levels "PDT000002821.3",..: 2 1 1 144 145 201 371 371 386 386 ...
    ##  $ Strand          : Factor w/ 2 levels "-","+": 1 2 2 1 1 1 2 2 2 2 ...
    ##  $ Element_symbol  : Factor w/ 1 level "blaSCO-1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Type            : Factor w/ 1 level "AMR": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Method          : Factor w/ 2 levels "ALLELEP","ALLELEX": 1 1 1 1 1 2 1 1 1 1 ...

### 2.3 Metadata Cleaning and Standardization

To facilitate downstream analysis, the Scientific_name column was parsed
to extract genus-level information, which was stored as a new factor
variable Genus. Additionally, all entries containing *Salmonella* were
standardized under a unified Species label: *Salmonella enterica*. Other
transformations are performed to enabled consistent grouping and
summarization across diverse isolate entries.

``` r
# Extract Genus
data$Genus <- word(data$Scientific_name, 1) # create a new column with just the genus name
data$Genus <- as.factor(data$Genus) # convert  the new Genus data to factor

# Standardize Salmonella names
unique(data$Scientific_name) # Check for unique scientific names
```

    ##  [1] Klebsiella pneumoniae      Salmonella Wien           
    ##  [3] Escherichia coli           Salmonella enterica       
    ##  [5] Pseudomonas aeruginosa     Klebsiella oxytoca        
    ##  [7] Salmonella  Concord        Klebsiella quasipneumoniae
    ##  [9] Providencia rettgeri       Proteus mirabilis         
    ## [11] Enterobacter hormaechei    Salmonella 4,[5],12:i:-   
    ## [13] Morganella morganii        Klebsiella variicola      
    ## [15] Klebsiella michiganensis   Citrobacter freundii      
    ## [17] Pseudomonas putida         Shigella flexneri         
    ## [19] Acinetobacter baumannii    Vibrio cholerae           
    ## [21] Salmonella Isangi         
    ## 21 Levels: Acinetobacter baumannii ... Vibrio cholerae

``` r
data$Species <- data$Scientific_name  # Create a species column
data$Species[grepl("Salmonella", data$Species)] <- "Salmonella enterica"  # Replace all Salmonella with Salmonella enterica


# # Fix column name for clarity
colnames(data)[which(names(data) == "Isolation._type")] <- "Isolation_type"


# Standardize and clean multiple columns in one loop
columns_to_clean <- c("Host", "Isolation_type", "Isolation_source")

for (col in columns_to_clean) {
  # Convert to character
  data[[col]] <- as.character(data[[col]])
  
  # Replace NA values
  data[[col]][is.na(data[[col]])] <- "unknown"
  
  # Trim whitespace and convert to lowercase
  data[[col]] <- trimws(tolower(data[[col]]))
  
  # Convert back to factor
  data[[col]] <- factor(data[[col]])
}


# View first six rows of cleaned data
head(data)
```

    ##         Scientific_name Isolation_type  Location    Isolation_source    Host
    ## 1 Klebsiella pneumoniae       clinical Argentina cerebrospinal fluid   human
    ## 2       Salmonella Wien       clinical    France               human   human
    ## 3       Salmonella Wien       clinical    France               human   human
    ## 4      Escherichia coli       clinical      <NA>             unknown unknown
    ## 5      Escherichia coli       clinical      <NA>             unknown unknown
    ## 6       Salmonella Wien       clinical    France             unknown   human
    ##   Collection_date        Isolate Strand Element_symbol Type  Method       Genus
    ## 1            1980 PDT000011670.2      -       blaSCO-1  AMR ALLELEP  Klebsiella
    ## 2            1988 PDT000002821.3      +       blaSCO-1  AMR ALLELEP  Salmonella
    ## 3            1988 PDT000002821.3      +       blaSCO-1  AMR ALLELEP  Salmonella
    ## 4            1988 PDT000570875.1      -       blaSCO-1  AMR ALLELEP Escherichia
    ## 5            1988 PDT000570876.1      -       blaSCO-1  AMR ALLELEP Escherichia
    ## 6            1988 PDT001024780.1      -       blaSCO-1  AMR ALLELEX  Salmonella
    ##                 Species
    ## 1 Klebsiella pneumoniae
    ## 2   Salmonella enterica
    ## 3   Salmonella enterica
    ## 4      Escherichia coli
    ## 5      Escherichia coli
    ## 6   Salmonella enterica

``` r
# Save the cleaned dataset
write.csv(data, "Data/blasco_cleaned_data.csv", row.names = FALSE)
```

## 3. Data Exploration

### 3.1 Data Overview

To understand the overall structure and scope of the cleaned dataset, a
summary statistics and table summarizing key variables such as isolate
count, species diversity, host range, geographic distribution, and
collection timeframe were generated.

``` r
# Check summary statistics
summary(data)
```

    ##                Scientific_name       Isolation_type         Location  
    ##  Klebsiella pneumoniae :529    agricultural :  1    South Africa:103  
    ##  Salmonella  Concord   : 87    clinical     :525    Spain       : 72  
    ##  Klebsiella variicola  : 36    environmental: 15    Ethiopia    : 54  
    ##  Salmonella enterica   : 25    unknown      :227    Kenya       : 33  
    ##  Pseudomonas aeruginosa: 18                         Nigeria     : 32  
    ##  Escherichia coli      : 12                         (Other)     :246  
    ##  (Other)               : 61                         NA's        :228  
    ##           Isolation_source      Host     Collection_date           Isolate   
    ##  unknown          :346     chicken:  8   Min.   :1980    PDT001665754.1:  3  
    ##  blood            :138     dog    :  3   1st Qu.:2013    PDT000002821.3:  2  
    ##  urine            : 61     horses : 10   Median :2016    PDT000047429.4:  2  
    ##  clinical material: 47     human  :504   Mean   :2015    PDT000047432.4:  2  
    ##  wound            : 26     pig    :  7   3rd Qu.:2019    PDT000136259.2:  2  
    ##  rectal           : 20     unknown:236   Max.   :2024    PDT000136260.2:  2  
    ##  (Other)          :130                   NA's   :240     (Other)       :755  
    ##  Strand   Element_symbol  Type         Method             Genus    
    ##  -:200   blaSCO-1:768    AMR:768   ALLELEP:701   Klebsiella  :575  
    ##  +:568                             ALLELEX: 67   Salmonella  :121  
    ##                                                  Pseudomonas : 22  
    ##                                                  Escherichia : 12  
    ##                                                  Vibrio      : 12  
    ##                                                  Enterobacter: 10  
    ##                                                  (Other)     : 16  
    ##                    Species   
    ##  Klebsiella pneumoniae :529  
    ##  Salmonella enterica   :121  
    ##  Klebsiella variicola  : 36  
    ##  Pseudomonas aeruginosa: 18  
    ##  Escherichia coli      : 12  
    ##  Vibrio cholerae       : 12  
    ##  (Other)               : 40

``` r
# Create summary table of key variables
summary_table <- data.frame(
  Variable = c("Number of Isolates", "Bacterial Species", "Host Types", "Countries", 
               "Isolation Sources", "Isolation Type", "Time Range"),
  Value = c(nrow(data), length(unique(data$Scientific_name)), 
            length(unique(data$Host[!is.na(data$Host)])), length(unique(data$Location[!is.na(data$Location)])),
            length(unique(data$Isolation_source[!is.na(data$Isolation_source)])),
             length(unique(data$Isolation_type[!is.na(data$Isolation_type)])),
            paste(min(data$Collection_date, na.rm=TRUE), "-", 
                  max(data$Collection_date, na.rm=TRUE)))
)

summary_table
```

    ##             Variable       Value
    ## 1 Number of Isolates         768
    ## 2  Bacterial Species          21
    ## 3         Host Types           6
    ## 4          Countries          46
    ## 5  Isolation Sources          40
    ## 6     Isolation Type           4
    ## 7         Time Range 1980 - 2024

``` r
# Save summary table
write.csv(summary_table, "Results/Statistics_Outputs/summary_table.csv", row.names = FALSE)
```

### 3.2 Host and Species Distribution

To examine the distribution of *blaSCO-1*-harboring bacterial species
across various host types and isolation type, bar plots were generated

``` r
# Plot 1: Host distribution of blaSCO-1-harboring bacterial species

# Create the plot with distinct colors for each Host
SpHst <- ggplot(data, aes(x = Species, fill = Host)) +
  geom_bar(stat = "count", position = "stack") +
  labs(
    title = " ",
    x = "Species",
    y = "Number of Isolates",
    fill = "Host"
  ) +
  scale_fill_manual(values = c(
    "chicken" = "firebrick",
    "dog" = "darkorange",
    "horses" = "green",
    "human" = "steelblue",
    "pig" = "gold",
    "unknown" = "gray"
  )) + # Assign distinct colors to each host category
  theme_minimal(base_size = 12) + # A clean theme with a readable font size
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11)) # Rotate x-axis labels if needed
SpHst
```

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# Plot 2: Isolation Type and Host distribution of blaSCO-1-harboring bacterial species

SpHsIt <- ggplot(data, aes(x = Species, fill = Host)) +
  geom_bar(position = "stack") +
  facet_wrap(~ Isolation_type) +
  labs(
    title = " ",
    x = "Species",
    y = "Number of Isolates",
    fill = "Host"
  ) +
  scale_fill_manual(values = c(
    "chicken" = "firebrick",
    "dog" = "darkorange",
    "horses" = "green",
    "human" = "steelblue",
    "pig" = "gold",
    "unknown" = "gray"
  )) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11))

SpHsIt
```

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
# Plot 3: Arrange both ggplot objects into a single figure
 
SpHst_SpHsIt  <- ggarrange(
  SpHst,  # First plot: SpHst
  SpHsIt,  # Second plot: SpHsIt
  labels = "auto",  # Automatically label the plots (A, B, C, etc.)
  nrow = 1,  # Arrange the plots in 3 rows
  ncol = 2,  # Arrange the plots in 1 column
  common.legend = TRUE,   # Share one legend
  legend = "right"        # Position of the shared legend
)


# Save manuscript ready combined figure

jpeg("Results/Figures/HostSpecies_distribution_plots.jpg", width = 16, height = 8, units = "in", res = 300)
print(SpHst_SpHsIt)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

#### 3.2.1 Interpretation

The majority of *blaSCO-1*-positive isolates were identified as
*Klebsiella pneumoniae*, particularly from human hosts. A smaller number
of isolates were recovered from animals (chickens, pigs, horses),
indicating that while *blaSCO-1* is present in both clinical and
zoonotic contexts, its burden is currently highest in humans.

In addition, clinical isolates were the dominant isolation type across
all species. Most animal-origin isolates came from environmental or
unknown sources, whereas human isolates were almost entirely clinical.
This suggests strong clinical surveillance but possible under-detection
in animals or the environment.

### 3.3 Geographic Distribution

To assess the global distribution of bacterial isolates carrying the
*blaSCO-1* gene, country-level isolate counts were mapped using
choropleth visualization. Country names were extracted from the Location
column of the dataset, standardized, and matched with base map data.

``` r
# Plot 4
# Prepare location data: remove missing entries and standardize country names
location_data <- data %>%
  filter(!is.na(Location)) %>%                 # Remove isolates with missing location data
  mutate(Country = Location) %>%               # Rename for clarity
  group_by(Country) %>%
  summarize(Count = n())                       # Count isolates per country


# Get map data for world countries
world_map <- map_data("world")

# Join isolate count data with geographic map data
world_data <- left_join(world_map, location_data, by = c("region" = "Country"))

# Replace missing count values with 0 to avoid plotting issues
world_data$Count[is.na(world_data$Count)] <- 0

# Plot choropleth map of isolate counts by country
Geo_dist <- ggplot(data = world_data) +
  geom_polygon(aes(x = long, y = lat, group = group, fill = Count),
               color = "grey", linewidth = 0.1) +
  scale_fill_gradientn(
  colors = c("lightblue", "yellow", "orange", "darkred"),  # Gradient from low to high isolate counts
  values = scales::rescale(c(0, 1, 5, max(world_data$Count, na.rm = TRUE))),  # Emphasizes lower-count variation
  limits = c(0, max(world_data$Count, na.rm = TRUE)),  # Ensures full data range is represented
  name = "No. of\nIsolates"  # Legend title with line break
  ) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  labs(title = " ") +
  coord_fixed(1.3)
Geo_dist
```

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Save manuscript ready figure
jpeg("Results/Figures/global_distribution_map.jpg", width = 12, height = 8, units = "in", res = 300)
print(Geo_dist)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

#### 3.3.1 Interpretation

Isolates harboring *blaSCO-1* were found in 46 countries, with South
Africa (103), Spain (72), and Ethiopia (54) showing the highest counts.
In addition, the wide geographic range reflects global dissemination
across both developed and developing nations.

### 3.4 Temporal Distribution

To explore the temporal dynamics of *blaSCO-1*-harboring isolates, the
overall yearly trends, species-specific distributions over time, and
their stratification by host were visualized. These plots help identify
patterns in the emergence and spread of this resistance gene across
different bacterial species and host types over the sampling period.

``` r
# Plot 5: Frequency of Temporal Distribution of blaSCO-1-Harboring Isolates 
Temp_data <- data %>%
  group_by(Collection_date) %>%
  summarise(count = n(), .groups = "drop")

Temp <- ggplot(Temp_data, aes(x = Collection_date, y = count)) +
  geom_line(color = "black", linewidth = 0.7) +
  geom_point(color = "black", size = 2) +
  theme_bw(base_size = 12) +
  labs(
    title = " ",
    x = "Year",
    y = "Number of Isolates"
  )

Temp
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# Plot 6 : Temporal Distribution of blaSCO-1-Harboring Isolates across Species
Temp_Sp <-ggplot(data, aes(x = Collection_date, fill = Species)) +
  geom_histogram(binwidth = 1, color = "black") +
  theme_bw(base_size = 12) +
  labs(
    title = " ",
    x = "Year",
    y = "Number of Isolates",
    fill = "Species"
  )
Temp_Sp
```

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
# Plot 7 : Temporal Distribution of blaSCO-1-Harboring Isolates across Species and Host
Temp_SpHst <- ggplot(data, aes(x = Collection_date, fill = Species)) +
  geom_histogram(binwidth = 1, color = "black") +
  facet_wrap(~ Host, ncol = 2) +
  theme_bw(base_size = 12) +
  labs(
    title = " ",
    x = "Year",
    y = "Number of Isolates",
    fill = "Species"
  )
Temp_SpHst
```

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
# Plot 8: Arrange both ggplot objects into a single figure
Temp_Sp_SpHst <- ggarrange(
  Temp,  # First plot: Temp
  Temp_Sp,  # Second plot: Temp_Sp
  Temp_SpHst,  # Third plot: Temp_SpHst
  labels = "auto",  # Automatically label the plots (A, B, C, etc.)
  nrow = 1,  # Arrange the plots in 3 rows
  ncol = 3,  # Arrange the plots in 1 column
  common.legend = TRUE,
  legend = TRUE  # Do not include a legend in the combined figure
)
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_line()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_line()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 240 rows containing non-finite outside the scale range (`stat_bin()`).
    ## Removed 240 rows containing non-finite outside the scale range (`stat_bin()`).

``` r
# Save manuscript ready combined figure
jpeg("Results/Figures/temporal_distribution_map.jpg", width = 12, height = 8, units = "in", res = 300)
print(Temp_Sp_SpHst)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

#### 3.4.1 Interpretation

The earliest detection of *blaSCO-1* dates back to 1980, but most
isolates were collected after 2013. There is a noticeable rise in
detections between 2016 and 2019, suggesting increased gene
dissemination.

*Klebsiella pneumoniae* was the first species found to harbor *blaSCO-1*
in 1980 and remained the most dominant throughout the dataset. Other
early detections included Salmonella enterica and Escherichia coli, both
of which also appeared in more recent years, along with additional
species. This trend reflects a gradual diversification of bacterial
hosts carrying the resistance gene.

Human-derived isolates span the entire 1980–2024 timeline, indicating
continuous clinical relevance. In contrast, animal-derived isolates
(from chickens, pigs) appear sporadically and are largely restricted to
post-2000 samples, while horse-derived isolates were detected primarily
before 2000. The presence of *blaSCO-1* in both human and animal
isolates over time suggests the potential for cross-species transmission

## 4. Statistical Analysis

### 4.1 Normality Test

To assess the distribution of isolate collection dates, a Shapiro-Wilk
test was first conducted to evaluate normality. A histogram with an
overlaid normal distribution curve was then generated to visually
inspect how well the data approximate a normal distribution.

``` r
# Perform the Shapiro-Wilk test for normality
shapiro.test(data$Collection_date)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  data$Collection_date
    ## W = 0.81821, p-value < 2.2e-16

``` r
# Create a histogram and overlay a normal distribution curve
Norm_plot <- ggplot(data, aes(x = Collection_date)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
  stat_function(fun = dnorm, args = list(mean = mean(data$Collection_date, na.rm = TRUE), sd = sd(data$Collection_date, na.rm = TRUE)), 
                color = "red", linewidth = 1) +
  labs(title = " ",
       x = "Collection Date",
       y = "Density") +
  theme_minimal()

Norm_plot
```

    ## Warning: The dot-dot notation (`..density..`) was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `after_stat(density)` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Save manuscript ready combined figure
jpeg("Results/Figures/normality_plot.jpg", width = 12, height = 8, units = "in", res = 300)
print(Norm_plot)
```

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

#### 4.1.1 Interpretation

The Shapiro-Wilk test indicated that the Collection_date variable is not
normally distributed (W = 0.818, p \< 2.2e-16). The histogram confirms
skewness, validating the use of non-parametric tests for temporal
comparisons.

### 4.2 Kruskal-Wallis & Wilcoxon Tests:

#### 4.2.1 Host vs. Collection Date

This section evaluates whether the year of collection for
*blaSCO-1*-harboring isolates differs significantly among host species.
A Kruskal-Wallis test was used to assess overall group differences,
followed by pairwise Wilcoxon tests to identify specific host pairs with
significant differences. The boxplot visualization highlights these
differences, offering a temporal perspective on the emergence and spread
of blaSCO-1 across animal and human hosts.

#### 4.2.2 Isolation Type vs. Collection Date

To assess if there are significant differences in the collection year
across different isolation types. This is useful if you’re interested in
understanding whether the timing of collection varies depending on the
source of isolation (e.g., clinical, non-clinical).

``` r
# Variables to compare with Collection date
vars <- c("Host", "Isolation_type")

# Custom colors for Host
host_colors <- c(
  "chicken" = "firebrick",
  "dog" = "darkorange",
  "horses" = "green",
  "human" = "steelblue",
  "pig" = "gold",
  "unknown" = "gray"
)

# Ensure Host factor levels match color order
data$Host <- factor(data$Host, levels = names(host_colors))

# Initialize plot list
plot_list <- list()


for (v in vars) {
  formula <- as.formula(paste("Collection_date ~", v))
  
  # Kruskal-Wallis test
  kruskal_result <- kruskal.test(formula, data = data)
  print(paste("Kruskal-Wallis test for", v, "p-value:", round(kruskal_result$p.value, 8)))

  # Initialize an empty list to store pairwise comparisons 
  my_comparisons <- list()

  # # If the variable is 'Host' and the Kruskal-Wallis test is significant, perform pairwise Wilcoxon test for Host
  if (v == "Host" && kruskal_result$p.value < 0.05) {
    pw <- pairwise.wilcox.test(data$Collection_date, data[[v]], p.adjust.method = "BH")
    print("Pairwise Wilcoxon adjusted p-values:")
    print(pw$p.value)

    # Extract significant comparisons
    p_values <- pw$p.value
    sig_pairs <- which(p_values < 0.05, arr.ind = TRUE)
    
# Create the list of comparisons properly
    if (nrow(sig_pairs) > 0) {
      for (i in 1:nrow(sig_pairs)) {
        my_comparisons[[i]] <- c(rownames(p_values)[sig_pairs[i, 1]], colnames(p_values)[sig_pairs[i, 2]])
      }
    }
  }

  # Get max collection year per group for proper y-label placement
  max_y <- data %>%
    group_by(.data[[v]]) %>%
    summarise(max_val = max(Collection_date, na.rm = TRUE)) %>%
    pull(max_val) %>%
    max(na.rm = TRUE)

  label_y_pos <- max_y + 10  # Adjust spacing

  # Build plot with jitter
  p <- ggplot(data, aes_string(x = v, y = "Collection_date", fill = v)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(" ") +
    stat_compare_means(method = "kruskal.test", label.y = label_y_pos)

  # Add custom fill and Wilcoxon for Host
  if (v == "Host") {
    p <- p + scale_fill_manual(values = host_colors)
    if (length(my_comparisons) > 0) {
      p <- p + stat_compare_means(
        comparisons = my_comparisons,
        method = "wilcox.test",
        p.adjust.method = "BH",
        label = "p.signif"
      )
    }
  }

  # Store
  plot_list[[v]] <- p
  print(p)
}
```

    ## [1] "Kruskal-Wallis test for Host p-value: 3.14e-06"

    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties

    ## [1] "Pairwise Wilcoxon adjusted p-values:"
    ##             chicken        dog       horses     human       pig
    ## dog     1.000000000         NA           NA        NA        NA
    ## horses  0.001618742 0.04711125           NA        NA        NA
    ## human   0.782207939 0.78295121 1.127903e-06        NA        NA
    ## pig     0.819549567 1.00000000 3.085576e-03 0.7822079        NA
    ## unknown 0.782207939 0.84798781 4.940640e-02 0.1329060 0.7829512

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_compare_means()`).

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_signif()`).

    ## Warning in wilcox.test.default(c(1989, 1989, 1991, 1991, 1995, 1995, 1995, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(1989, 1989, 1991, 1991, 1995, 1995, 1995, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(2014, 2014, 2014, 2016, 2016, 2016, 2016:
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(1988, 1988, 1992, 2011, 2011, 2013, 2014, :
    ## cannot compute exact p-value with ties

    ## Warning: Removed 240 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

    ## [1] "Kruskal-Wallis test for Isolation_type p-value: 0.77156424"

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_compare_means()`).

    ## Warning: Removed 240 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
# Arrange both plots
Col <- ggarrange(
  plot_list[["Host"]],
  plot_list[["Isolation_type"]],
  labels = c("a", "b"),
  ncol = 2,
  common.legend = FALSE
)
```

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_compare_means()`).

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_signif()`).

    ## Warning in wilcox.test.default(c(1989, 1989, 1991, 1991, 1995, 1995, 1995, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(1989, 1989, 1991, 1991, 1995, 1995, 1995, :
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(2014, 2014, 2014, 2016, 2016, 2016, 2016:
    ## cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(1988, 1988, 1992, 2011, 2011, 2013, 2014, :
    ## cannot compute exact p-value with ties

    ## Warning: Removed 240 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Warning: Removed 240 rows containing non-finite outside the scale range
    ## (`stat_compare_means()`).

    ## Warning: Removed 240 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
# Save manuscript ready combined figure
jpeg("Results/Figures/collection_date_significance.jpg", width = 12, height = 8, units = "in", res = 300)
print(Col)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# Save wilcoxon pairwise result for Host Comparison
write.csv(pw$p.value, "Results/Statistics_Outputs/pairwise_wilcox_Host.csv")
```

#### 4.2.3 Interpretation

The Kruskal-Wallis test showed a statistically significant difference in
the collection year across host types (p = 3.14e-06). Specifically,
isolates from horses were significantly older compared to those from
pigs, humans, chicken and unknown hosts, suggesting either earlier
emergence of this group.

There was no significant difference in collection year between different
isolation types (p = 0.77). This suggests that the temporal distribution
of blaSCO-1 is consistent regardless of whether the isolate was
clinical, environmental, or from agricultural settings.

### 4.3 Fisher’s Test:

#### 4.3.1 Isolation source vs. Host

To evaluate whether there is a significant association between the
source of isolation and the host species of *blaSCO-1*-positive
isolates. It helps determine if certain hosts are more likely to be
associated with specific isolation sources (e.g., human samples from
wounds or animals from feces).

``` r
# Create a contingency table for Isolation Source vs Host
table_isolation_host <- table(data$Isolation_source, data$Host)
head(table_isolation_host)
```

    ##                     
    ##                      chicken dog horses human pig unknown
    ##   bile                     0   0      0     2   0       0
    ##   biopsy                   0   0      0     2   0       0
    ##   blood                    0   0      0   138   0       0
    ##   bronchial aspirate       0   0      0     1   0       0
    ##   burn                     0   0      0     1   0       1
    ##   catheter                 0   0      0     3   0       0

``` r
# Fisher's Exact test will be used because many isolate source count were zero(0) for some hosts

# Perform Fisher's Exact Test using simulation
fisher_result_isolation_host_sim <- fisher.test(table_isolation_host, simulate.p.value = TRUE, B = 10000)

# Print the result
print(fisher_result_isolation_host_sim)
```

    ## 
    ##  Fisher's Exact Test for Count Data with simulated p-value (based on
    ##  10000 replicates)
    ## 
    ## data:  table_isolation_host
    ## p-value = 9.999e-05
    ## alternative hypothesis: two.sided

``` r
Is_Hst <- ggplot(data, aes(x = Isolation_source, y = Host, color = Host)) +
  geom_jitter(width = 0.3, height = 0.2, alpha = 0.8, size = 1) +
  labs(title = " ",
       x = "Isolation Source",
       y = "Host") +
  theme_minimal() +
  scale_color_manual(values = host_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip() # flip axis
Is_Hst
```

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# Save manuscript ready figure
jpeg("Results/Figures/Host_isolate_source_association.jpg", width = 12, height = 8, units = "in", res = 300)
print(Is_Hst)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

#### 4.3.2 Interpretation

Fisher’s Exact Test revealed a significant association between isolation
source and host species (simulated p-value = 9.999e-05). This suggests
that specific hosts are more likely to be associated with particular
types of samples. For example, human isolates were primarily recovered
from clinical sources such as urine, blood, wound, while isolates from
animals (e.g., chickens, pigs) were more frequently linked to feces. The
concentration of human-derived isolates in clinical materials reinforces
the gene’s clinical relevance, whereas the sparse distribution and
concentration in “unknown” source in animal hosts indicates either
limited surveillance or poor documentation in non-human hosts

### 4.4 Chi-Square Test:

#### 4.4.1 Isolation Type vs. Presence of Gene

To test whether the presence of the *blaSCO-1* gene is associated with
being isolated from clinical vs. non-clinical sources. This helps
understand whether there is a relationship between the gene’s presence
and the origin of the sample.

``` r
# Create a contingency table for clinical vs. non-clinical isolates
table_clinical_nonclinical <- table(data$Isolation_type, data$Element_symbol) 
print(table_clinical_nonclinical)
```

    ##                
    ##                 blaSCO-1
    ##   agricultural         1
    ##   clinical           525
    ##   environmental       15
    ##   unknown            227

``` r
# Perform Chi-square test
chi_square_result <- chisq.test(table_clinical_nonclinical)
print(chi_square_result)
```

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  table_clinical_nonclinical
    ## X-squared = 937.1, df = 3, p-value < 2.2e-16

#### 4.4.2 Interpretation

The chi-square test indicated a highly significant association between
the presence of *blaSCO-1* and clinical origin (p \< 2.2e-16). Over 68%
of isolates were clinical, reinforcing the gene’s strong association
with healthcare settings.

### 4.5 Generalized Linear Model:

#### 4.5.1 Species vs. Host Interaction

To investigate how bacterial species, host type, and their interaction
influence the year of detection for *blaSCO-1*-harboring isolates, a
generalized linear model (GLM) was fitted with Collection_date as the
response variable. This approach allows us to assess whether certain
species or host groups are associated with earlier or more recent
detections, and whether combinations of the two explain additional
variation even though our data is not normally distributed

``` r
# Fit the generalized linear model
glm_model <- glm(Collection_date ~ Species * Host, data = data)
summary(glm_model)
```

    ## 
    ## Call:
    ## glm(formula = Collection_date ~ Species * Host, data = data)
    ## 
    ## Coefficients: (60 not defined because of singularities)
    ##                                                 Estimate Std. Error t value
    ## (Intercept)                                    2.032e+03  6.984e+00 290.910
    ## SpeciesEnterobacter hormaechei                -1.982e+01  8.437e+00  -2.349
    ## SpeciesEscherichia coli                       -3.432e+01  7.394e+00  -4.641
    ## SpeciesKlebsiella michiganensis               -3.253e-12  6.635e+00   0.000
    ## SpeciesKlebsiella oxytoca                     -1.600e+01  6.635e+00  -2.411
    ## SpeciesKlebsiella pneumoniae                  -1.682e+01  6.784e+00  -2.479
    ## SpeciesKlebsiella quasipneumoniae             -1.043e+01  5.016e+00  -2.079
    ## SpeciesKlebsiella variicola                   -3.459e-12  4.757e+00   0.000
    ## SpeciesMorganella morganii                    -2.500e+00  5.746e+00  -0.435
    ## SpeciesProteus mirabilis                      -1.250e-01  4.977e+00  -0.025
    ## SpeciesProvidencia rettgeri                   -6.000e+00  6.635e+00  -0.904
    ## SpeciesPseudomonas aeruginosa                 -4.182e+01  8.437e+00  -4.956
    ## SpeciesPseudomonas putida                      5.000e+00  5.246e+00   0.953
    ## SpeciesSalmonella enterica                    -1.682e+01  4.901e+00  -3.432
    ## Hostdog                                       -7.255e-13  3.176e+00   0.000
    ## Hosthorses                                    -2.140e+01  5.193e+00  -4.121
    ## Hosthuman                                     -1.282e+01  5.174e+00  -2.478
    ## Hostpig                                        1.000e+00  3.709e+00   0.270
    ## Hostunknown                                    2.000e+00  2.428e+00   0.824
    ## SpeciesEnterobacter hormaechei:Hostdog                NA         NA      NA
    ## SpeciesEscherichia coli:Hostdog                       NA         NA      NA
    ## SpeciesKlebsiella michiganensis:Hostdog               NA         NA      NA
    ## SpeciesKlebsiella oxytoca:Hostdog                     NA         NA      NA
    ## SpeciesKlebsiella pneumoniae:Hostdog                  NA         NA      NA
    ## SpeciesKlebsiella quasipneumoniae:Hostdog             NA         NA      NA
    ## SpeciesKlebsiella variicola:Hostdog                   NA         NA      NA
    ## SpeciesMorganella morganii:Hostdog                    NA         NA      NA
    ## SpeciesProteus mirabilis:Hostdog                      NA         NA      NA
    ## SpeciesProvidencia rettgeri:Hostdog                   NA         NA      NA
    ## SpeciesPseudomonas aeruginosa:Hostdog                 NA         NA      NA
    ## SpeciesPseudomonas putida:Hostdog                     NA         NA      NA
    ## SpeciesSalmonella enterica:Hostdog                    NA         NA      NA
    ## SpeciesEnterobacter hormaechei:Hosthorses             NA         NA      NA
    ## SpeciesEscherichia coli:Hosthorses                    NA         NA      NA
    ## SpeciesKlebsiella michiganensis:Hosthorses            NA         NA      NA
    ## SpeciesKlebsiella oxytoca:Hosthorses                  NA         NA      NA
    ## SpeciesKlebsiella pneumoniae:Hosthorses               NA         NA      NA
    ## SpeciesKlebsiella quasipneumoniae:Hosthorses          NA         NA      NA
    ## SpeciesKlebsiella variicola:Hosthorses                NA         NA      NA
    ## SpeciesMorganella morganii:Hosthorses                 NA         NA      NA
    ## SpeciesProteus mirabilis:Hosthorses                   NA         NA      NA
    ## SpeciesProvidencia rettgeri:Hosthorses                NA         NA      NA
    ## SpeciesPseudomonas aeruginosa:Hosthorses              NA         NA      NA
    ## SpeciesPseudomonas putida:Hosthorses                  NA         NA      NA
    ## SpeciesSalmonella enterica:Hosthorses                 NA         NA      NA
    ## SpeciesEnterobacter hormaechei:Hosthuman       1.915e+01  7.185e+00   2.666
    ## SpeciesEscherichia coli:Hosthuman              2.732e+01  7.394e+00   3.694
    ## SpeciesKlebsiella michiganensis:Hosthuman             NA         NA      NA
    ## SpeciesKlebsiella oxytoca:Hosthuman                   NA         NA      NA
    ## SpeciesKlebsiella pneumoniae:Hosthuman         1.282e+01  4.906e+00   2.614
    ## SpeciesKlebsiella quasipneumoniae:Hosthuman           NA         NA      NA
    ## SpeciesKlebsiella variicola:Hosthuman                 NA         NA      NA
    ## SpeciesMorganella morganii:Hosthuman                  NA         NA      NA
    ## SpeciesProteus mirabilis:Hosthuman                    NA         NA      NA
    ## SpeciesProvidencia rettgeri:Hosthuman                 NA         NA      NA
    ## SpeciesPseudomonas aeruginosa:Hosthuman        4.176e+01  7.104e+00   5.878
    ## SpeciesPseudomonas putida:Hosthuman                   NA         NA      NA
    ## SpeciesSalmonella enterica:Hosthuman                  NA         NA      NA
    ## SpeciesEnterobacter hormaechei:Hostpig                NA         NA      NA
    ## SpeciesEscherichia coli:Hostpig                1.550e+01  5.196e+00   2.983
    ## SpeciesKlebsiella michiganensis:Hostpig               NA         NA      NA
    ## SpeciesKlebsiella oxytoca:Hostpig                     NA         NA      NA
    ## SpeciesKlebsiella pneumoniae:Hostpig                  NA         NA      NA
    ## SpeciesKlebsiella quasipneumoniae:Hostpig             NA         NA      NA
    ## SpeciesKlebsiella variicola:Hostpig                   NA         NA      NA
    ## SpeciesMorganella morganii:Hostpig                    NA         NA      NA
    ## SpeciesProteus mirabilis:Hostpig                      NA         NA      NA
    ## SpeciesProvidencia rettgeri:Hostpig                   NA         NA      NA
    ## SpeciesPseudomonas aeruginosa:Hostpig                 NA         NA      NA
    ## SpeciesPseudomonas putida:Hostpig                     NA         NA      NA
    ## SpeciesSalmonella enterica:Hostpig                    NA         NA      NA
    ## SpeciesEnterobacter hormaechei:Hostunknown            NA         NA      NA
    ## SpeciesEscherichia coli:Hostunknown                   NA         NA      NA
    ## SpeciesKlebsiella michiganensis:Hostunknown           NA         NA      NA
    ## SpeciesKlebsiella oxytoca:Hostunknown                 NA         NA      NA
    ## SpeciesKlebsiella pneumoniae:Hostunknown              NA         NA      NA
    ## SpeciesKlebsiella quasipneumoniae:Hostunknown         NA         NA      NA
    ## SpeciesKlebsiella variicola:Hostunknown               NA         NA      NA
    ## SpeciesMorganella morganii:Hostunknown                NA         NA      NA
    ## SpeciesProteus mirabilis:Hostunknown                  NA         NA      NA
    ## SpeciesProvidencia rettgeri:Hostunknown               NA         NA      NA
    ## SpeciesPseudomonas aeruginosa:Hostunknown             NA         NA      NA
    ## SpeciesPseudomonas putida:Hostunknown                 NA         NA      NA
    ## SpeciesSalmonella enterica:Hostunknown                NA         NA      NA
    ##                                               Pr(>|t|)    
    ## (Intercept)                                    < 2e-16 ***
    ## SpeciesEnterobacter hormaechei                0.019215 *  
    ## SpeciesEscherichia coli                       4.43e-06 ***
    ## SpeciesKlebsiella michiganensis               1.000000    
    ## SpeciesKlebsiella oxytoca                     0.016251 *  
    ## SpeciesKlebsiella pneumoniae                  0.013504 *  
    ## SpeciesKlebsiella quasipneumoniae             0.038112 *  
    ## SpeciesKlebsiella variicola                   1.000000    
    ## SpeciesMorganella morganii                    0.663708    
    ## SpeciesProteus mirabilis                      0.979971    
    ## SpeciesProvidencia rettgeri                   0.366297    
    ## SpeciesPseudomonas aeruginosa                 9.83e-07 ***
    ## SpeciesPseudomonas putida                     0.340966    
    ## SpeciesSalmonella enterica                    0.000649 ***
    ## Hostdog                                       1.000000    
    ## Hosthorses                                    4.41e-05 ***
    ## Hosthuman                                     0.013555 *  
    ## Hostpig                                       0.787583    
    ## Hostunknown                                   0.410543    
    ## SpeciesEnterobacter hormaechei:Hostdog              NA    
    ## SpeciesEscherichia coli:Hostdog                     NA    
    ## SpeciesKlebsiella michiganensis:Hostdog             NA    
    ## SpeciesKlebsiella oxytoca:Hostdog                   NA    
    ## SpeciesKlebsiella pneumoniae:Hostdog                NA    
    ## SpeciesKlebsiella quasipneumoniae:Hostdog           NA    
    ## SpeciesKlebsiella variicola:Hostdog                 NA    
    ## SpeciesMorganella morganii:Hostdog                  NA    
    ## SpeciesProteus mirabilis:Hostdog                    NA    
    ## SpeciesProvidencia rettgeri:Hostdog                 NA    
    ## SpeciesPseudomonas aeruginosa:Hostdog               NA    
    ## SpeciesPseudomonas putida:Hostdog                   NA    
    ## SpeciesSalmonella enterica:Hostdog                  NA    
    ## SpeciesEnterobacter hormaechei:Hosthorses           NA    
    ## SpeciesEscherichia coli:Hosthorses                  NA    
    ## SpeciesKlebsiella michiganensis:Hosthorses          NA    
    ## SpeciesKlebsiella oxytoca:Hosthorses                NA    
    ## SpeciesKlebsiella pneumoniae:Hosthorses             NA    
    ## SpeciesKlebsiella quasipneumoniae:Hosthorses        NA    
    ## SpeciesKlebsiella variicola:Hosthorses              NA    
    ## SpeciesMorganella morganii:Hosthorses               NA    
    ## SpeciesProteus mirabilis:Hosthorses                 NA    
    ## SpeciesProvidencia rettgeri:Hosthorses              NA    
    ## SpeciesPseudomonas aeruginosa:Hosthorses            NA    
    ## SpeciesPseudomonas putida:Hosthorses                NA    
    ## SpeciesSalmonella enterica:Hosthorses               NA    
    ## SpeciesEnterobacter hormaechei:Hosthuman      0.007932 ** 
    ## SpeciesEscherichia coli:Hosthuman             0.000244 ***
    ## SpeciesKlebsiella michiganensis:Hosthuman           NA    
    ## SpeciesKlebsiella oxytoca:Hosthuman                 NA    
    ## SpeciesKlebsiella pneumoniae:Hosthuman        0.009226 ** 
    ## SpeciesKlebsiella quasipneumoniae:Hosthuman         NA    
    ## SpeciesKlebsiella variicola:Hosthuman               NA    
    ## SpeciesMorganella morganii:Hosthuman                NA    
    ## SpeciesProteus mirabilis:Hosthuman                  NA    
    ## SpeciesProvidencia rettgeri:Hosthuman               NA    
    ## SpeciesPseudomonas aeruginosa:Hosthuman       7.55e-09 ***
    ## SpeciesPseudomonas putida:Hosthuman                 NA    
    ## SpeciesSalmonella enterica:Hosthuman                NA    
    ## SpeciesEnterobacter hormaechei:Hostpig              NA    
    ## SpeciesEscherichia coli:Hostpig               0.002989 ** 
    ## SpeciesKlebsiella michiganensis:Hostpig             NA    
    ## SpeciesKlebsiella oxytoca:Hostpig                   NA    
    ## SpeciesKlebsiella pneumoniae:Hostpig                NA    
    ## SpeciesKlebsiella quasipneumoniae:Hostpig           NA    
    ## SpeciesKlebsiella variicola:Hostpig                 NA    
    ## SpeciesMorganella morganii:Hostpig                  NA    
    ## SpeciesProteus mirabilis:Hostpig                    NA    
    ## SpeciesProvidencia rettgeri:Hostpig                 NA    
    ## SpeciesPseudomonas aeruginosa:Hostpig               NA    
    ## SpeciesPseudomonas putida:Hostpig                   NA    
    ## SpeciesSalmonella enterica:Hostpig                  NA    
    ## SpeciesEnterobacter hormaechei:Hostunknown          NA    
    ## SpeciesEscherichia coli:Hostunknown                 NA    
    ## SpeciesKlebsiella michiganensis:Hostunknown         NA    
    ## SpeciesKlebsiella oxytoca:Hostunknown               NA    
    ## SpeciesKlebsiella pneumoniae:Hostunknown            NA    
    ## SpeciesKlebsiella quasipneumoniae:Hostunknown       NA    
    ## SpeciesKlebsiella variicola:Hostunknown             NA    
    ## SpeciesMorganella morganii:Hostunknown              NA    
    ## SpeciesProteus mirabilis:Hostunknown                NA    
    ## SpeciesProvidencia rettgeri:Hostunknown             NA    
    ## SpeciesPseudomonas aeruginosa:Hostunknown           NA    
    ## SpeciesPseudomonas putida:Hostunknown               NA    
    ## SpeciesSalmonella enterica:Hostunknown              NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 22.014)
    ## 
    ##     Null deviance: 20822  on 527  degrees of freedom
    ## Residual deviance: 11095  on 504  degrees of freedom
    ##   (240 observations deleted due to missingness)
    ## AIC: 3156.2
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
# Estimate EMMeans
emm <- emmeans(glm_model, ~ Species * Host)

# Convert to data frame and remove rows with missing emmeans (non-estimable)
emm_df <- as.data.frame(emm)
emm_df <- na.omit(emm_df)

# Plot only the valid combinations
library(ggplot2)
em_plot <-ggplot(emm_df, aes(x = interaction(Species, Host), y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = " ",
    x = "Species–Host Combination",
    y = "Estimated Mean Collection Year"
  )
em_plot
```

![](Blasco-1-Analysis_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# Save manuscript ready figure
jpeg("Results/Figures/Host_species_interaction.jpg", width = 12, height = 8, units = "in", res = 300)
print(em_plot)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# Save Estimated Marginal Means of species and host interaction
write.csv(emm_df, "Results/Statistics_Outputs/emmeans_species_host.csv", row.names = FALSE)
```

#### 4.5.2 Interpretation

A generalized linear model (GLM) was used to assess how Species, Host,
and their interaction influenced the year of isolate collection. The
model explained a moderate portion of the variation in the data
(adjusted R² ≈ 0.44; AIC = 3156.2). Significant predictors included both
main effects (species and host) and several species–host combinations,
although 60 interaction terms were not estimable were not estimable due
to singularities or limited temporal variation within these groups

Key findings:

*Klebsiella pneumoniae*, *Escherichia coli*, *Pseudomonas aeruginosa*,
and *Salmonella enterica* were all significantly associated with earlier
or later collection dates compared to the baseline.

Host species such as horses (p = 4.41e-05) and humans (p = 0.0136)
showed significant effects, with horse-derived isolates collected
significantly earlier.

Several meaningful interactions were identified:

*E. coli* in humans and pigs was associated with more recent collection
dates (p = 0.00024 and p = 0.0030, respectively).

*P. aeruginosa* in humans also showed a strong effect (p = 7.55e-09),
consistent with clinical emergence in later years.

*K. pneumoniae* in humans was significantly associated with earlier
detections (p = 0.0092), aligning with its appearance in 1980.

These results reinforce the temporal stratification of blaSCO-1 by both
species and host, with implications for its historical emergence and
current distribution in clinical versus animal reservoirs.

## 5. Discussion and Conclusion

This study presents the most extensive exploration to date of the
distribution, host range, and temporal emergence of bacterial isolates
harboring the *blaSCO-1* resistance gene, using metadata from the NCBI
Pathogen Detection Database. Despite the our initial report of the
gene’s presence in *Salmonella* isolated decades ago, our findings
demonstrate that *blaSCO-1* has since disseminated across multiple
genera and continents, suggesting that it is more widespread than
previously recognized.

A striking pattern emerged from the temporal analysis, showing that
while *blaSCO-1* was first detected in 1980, the majority of isolates
have been reported since 2013, with a clear surge between 2016 and 2019.
This likely reflects both increasing gene dissemination and improved
detection methods in recent years. *Klebsiella pneumoniae* was
identified as the earliest and most dominant host, though other species
such as *Escherichia coli*, *Pseudomonas aeruginosa*, and *Salmonella
enterica* were also significantly represented—particularly in later
years. This supports a hypothesis of progressive horizontal gene
transfer and ecological expansion.

The gene was strongly associated with human clinical isolates, as
confirmed by the chi-square test (p \< 2.2e-16), which revealed that
over 68% of *blaSCO-1*-harboring isolates originated from clinical
settings. This suggests that *blaSCO-1* is primarily a threat in
healthcare environments but may also be underreported in animals and
environmental sources. The Fisher’s Exact Test further confirmed a
significant relationship between host type and isolation source, with
human isolates linked to specific clinical materials (e.g., blood,
urine), while animal isolates were mostly associated with feces or
labeled as “unknown.” These patterns may reflect both real biological
distributions and limitations in surveillance outside human health.

Host–species modeling revealed that horses were disproportionately
associated with earlier detections, while *P. aeruginosa* and *E. coli*
in humans and pigs were tied to more recent isolates. However,
limitations in statistical inference were noted for many species–host
pairs due to data sparsity or singularities in the linear model,
particularly in under-sampled groups such as dogs or agricultural
settings.

Together, these results emphasize the need for enhanced One Health
surveillance especially because *blaSCO-1* may be circulating silently.
The zoonotic potential of this resistance gene cannot be ignored, given
its presence across both human and animal hosts. Moreover, its
association with extended-spectrum β-lactam resistance raises concerns
for therapeutic efficacy in both sectors.

In conclusion, *blaSCO-1* represents an underrecognized resistance gene
with broad host range and global distribution. Its emergence in clinical
settings is has not been well documented, and gaps remain in our
understanding of its ecological reservoirs and transmission dynamics.
Further genomic studies and environmental surveillance are warranted to
fully characterize its risk and guide targeted intervention strategies.

## 6. Reference

[Iduu, N. V., Raiford, D., Conley, A., Scaria, J., Nelson, J., Ruesch,
L., Price, S., Yue, M., Gong, J., Wei, L., & Wang, C. (2024). A
Retrospective Analysis of Salmonella Isolates across 11 Animal Species
(1982-1999) Led to the First Identification of Chromosomally Encoded
blaSCO-1 in the USA. Microorganisms, 12(3),
528](%22https://doi.org/10.3390/microorganisms12030528%22)

## 7. Files Availability

[Link to Github
repository](https://github.com/NVI0001/blasco-1-analysis)
