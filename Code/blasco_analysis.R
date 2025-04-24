## ----setup, include=TRUE----------------------------------------------------------------------

# Load required libraries
library(stringr)      # For string manipulation functions (e.g., pattern matching, replacement, extraction)
library(dplyr)        # For data manipulation (e.g., filter, mutate, select, summarise, etc.)
library(tidyverse)    # Loads a collection of tidy data tools (includes dplyr, ggplot2, tidyr, readr, etc.)
library(ggplot2)      # For creating data visualizations and plots
library(forcats)      # For handling and reordering categorical/factor variables
library(maps)         # For drawing geographical maps (includes country/state boundaries)
library(lme4)         # For fitting linear and generalized linear mixed-effects models
library(ggpubr)       # For publication-ready plots, adds enhancements to ggplot2
library(ggrepel)      # For better label placement in ggplot2 (avoids label overlap)
library(emmeans)      # For estimating marginal means (EMMs) aka least-squares means from model outputs
library(multcomp)     # For multiple comparison procedures and simultaneous inference

install.packages("tinytex", repos = "https://cran.rstudio.com")  # Installs 'tinytex' for PDF rendering via LaTeX


#load data 
data <- read.csv("Data/blasco_rawdata.csv", na.strings = c("", "na"), stringsAsFactors = TRUE)

str(data) # View data structure



## ---------------------------------------------------------------------------------------------

# Extract Genus
data$Genus <- word(data$Scientific_name, 1) # create a new column with just the genus name
data$Genus <- as.factor(data$Genus) # convert  the new Genus data to factor

# Standardize Salmonella names
unique(data$Scientific_name) # Check for unique scientific names
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

# Save the cleaned dataset
write.csv(data, "Data/blasco_cleaned_data.csv", row.names = FALSE)


## ---------------------------------------------------------------------------------------------
# Check summary statistics
summary(data)

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

# Save summary table
write.csv(summary_table, "Data/summary_table.csv", row.names = FALSE)




## ---------------------------------------------------------------------------------------------

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




## ---------------------------------------------------------------------------------------------
# Plot 4: Geographic distribution of Isolates
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

# Save manuscript ready figure
jpeg("Results/Figures/global_distribution_map.jpg", width = 12, height = 8, units = "in", res = 300)
print(Geo_dist)
dev.off()


## ---------------------------------------------------------------------------------------------

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


# Save manuscript ready combined figure
jpeg("Results/Figures/temporal_distribution_map.jpg", width = 12, height = 8, units = "in", res = 300)
print(Temp_Sp_SpHst)
dev.off()


## ---------------------------------------------------------------------------------------------

# Perform the Shapiro-Wilk test for normality
shapiro.test(data$Collection_date)

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

# Save manuscript ready combined figure
jpeg("Results/Figures/normality_plot.jpg", width = 12, height = 8, units = "in", res = 300)
print(Norm_plot)
dev.off()


## ---------------------------------------------------------------------------------------------

# Kruskal-Wallis and Wilcoxon test: Variables to compare with Collection date
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

# Arrange both plots
Col <- ggarrange(
  plot_list[["Host"]],
  plot_list[["Isolation_type"]],
  labels = c("a", "b"),
  ncol = 2,
  common.legend = FALSE
)


# Save manuscript ready combined figure
jpeg("Results/Figures/collection_date_significance.jpg", width = 12, height = 8, units = "in", res = 300)
print(Col)
dev.off()


## ---------------------------------------------------------------------------------------------

# Create a contingency table for Isolation Source vs Host Association
table_isolation_host <- table(data$Isolation_source, data$Host)
head(table_isolation_host)

# Fisher's Exact test will be used because many isolate source count were zero(0) for some hosts

# Perform Fisher's Exact Test using simulation
fisher_result_isolation_host_sim <- fisher.test(table_isolation_host, simulate.p.value = TRUE, B = 10000)

# Print the result
print(fisher_result_isolation_host_sim)


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

# Save manuscript ready figure
jpeg("Results/Figures/Host_isolate_source_association.jpg", width = 12, height = 8, units = "in", res = 300)
print(Is_Hst)
dev.off()


## ---------------------------------------------------------------------------------------------

# Create a contingency table for clinical vs. non-clinical isolates
table_clinical_nonclinical <- table(data$Isolation_type, data$Element_symbol) 
print(table_clinical_nonclinical)


# Perform Chi-square test
chi_square_result <- chisq.test(table_clinical_nonclinical)
print(chi_square_result)



## ---------------------------------------------------------------------------------------------

# Fit the generalized linear model to see host vs species interaction over time
glm_model <- glm(Collection_date ~ Species * Host, data = data)
summary(glm_model)


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

# Save manuscript ready figure
jpeg("Results/Figures/Host_species_interaction.jpg", width = 12, height = 8, units = "in", res = 300)
print(em_plot)
dev.off()

