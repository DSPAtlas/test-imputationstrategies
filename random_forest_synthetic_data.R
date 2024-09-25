# Load required libraries
library(ggplot2)
library(ggpubr)
library(missForest)
library(dplyr)
library(tidyr)
library(protti)

# Function to generate synthetic data and process it
generate_and_process_data <- function(dropout_curve_inflection_value, dropout_curve_sd_value) {
  
  synthetic_data <- create_synthetic_data(
    n_proteins = 100,
    frac_change = 0.1,
    n_replicates = 3,
    n_conditions = 3,
    method = "effect_random", 
    dropout_curve_inflection = dropout_curve_inflection_value,
    dropout_curve_sd = dropout_curve_sd_value
  )
  
  # Processing the synthetic data
  synthetic_data_wide <- synthetic_data %>%
    tidyr::pivot_wider(
      id_cols=sample,
      names_from = peptide, 
      values_from = peptide_intensity_missing) %>%
    dplyr::select(!sample)
  
  synthetic_data_wide <- data.frame(lapply(synthetic_data_wide, function(x) as.numeric(as.character(x))))
  
  synthetic_data_imputed_rf <- missForest(synthetic_data_wide, verbose  = TRUE)
  synthetic_data_imputed <- synthetic_data_imputed_rf$ximp
  synthetic_data_imputed["sample"] <- synthetic_data %>%
    tidyr::pivot_wider(
      id_cols=sample,
      names_from = peptide, 
      values_from = peptide_intensity_missing) %>%
    dplyr::pull(sample)
  
  synthetic_data_imputed <- synthetic_data_imputed %>%
    as.data.frame() %>%                      
    tidyr::pivot_longer(
      cols = -sample,                    
      names_to = "peptide",           
      values_to = "peptide_intensity_imputed"
    )
  print(synthetic_data_imputed %>% nrow())
  
  synthetic_data_imputed <- synthetic_data_imputed %>%
    left_join(
      synthetic_data %>% 
        dplyr::select(sample, condition, peptide, protein, peptide_intensity,peptide_intensity_missing) %>% 
        distinct(),  
      by = c("sample", "peptide")
    ) %>%
    dplyr::filter(is.na(peptide_intensity_missing))
  print(synthetic_data_imputed %>% nrow())
  return(synthetic_data_imputed)
}

# Define dropout curve inflection and SD values
dropout_curve_inflection_values <-  c(10, 14, 20)
dropout_curve_sd_values <- c(-0.8,-1.2, -2)

# Create all combinations of inflection and sd values
combinations <- expand.grid(
  inflection = dropout_curve_inflection_values,
  sd_value = dropout_curve_sd_values
)

# Number of datasets based on combinations
num_datasets <- 30
datasets <- list()

# Loop through all combinations of dropout_curve_inflection_values and dropout_curve_sd_values
for (i in 1:9) {
  inflection_value <- combinations$inflection[i]
  sd_value <- combinations$sd_value[i]
  
  datasets[[i]] <- generate_and_process_data(
    dropout_curve_inflection_value = inflection_value,
    dropout_curve_sd_value = sd_value
  )
}

# Plot each dataset in a scatterplot
plot_list <- list()

for (i in 1:9) {
  # Calculate Pearson correlation and p-value
  cor_test <- cor.test(datasets[[i]]$peptide_intensity, datasets[[i]]$peptide_intensity_imputed, method = "pearson", use = "complete.obs")
  r_value <- round(cor_test$estimate, 2)
  p_value <- cor_test$p.value
  
  # Create plot
  plot <- datasets[[i]] %>%
    dplyr::filter(is.na(peptide_intensity_missing)) %>%
    ggplot(aes(x = peptide_intensity, y = peptide_intensity_imputed)) +
    geom_point(color = "blue") +              
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    labs(
      title = paste("Scatterplot for Dataset", i),
      subtitle = paste("Pearson r =", r_value, "| p-value =", p_value),  # Add r and p-value in the subtitle
      x = "Original Peptide Intensity",
      y = "Imputed Peptide Intensity"
    ) +
    xlim(min(datasets[[i]]$peptide_intensity, na.rm = TRUE), 
         max(datasets[[i]]$peptide_intensity, na.rm = TRUE)) +
    theme_minimal()
  
  plot_list[[i]] <- plot
}

# Display all plots together
ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)


