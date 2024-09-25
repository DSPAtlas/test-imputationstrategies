library(protti)
library(dplyr)
library(magrittr)
library(ggplot2)

# load data
utils::data("rapamycin_dose_response")

data_normalised <- rapamycin_dose_response %>%
  filter(eg_is_decoy == FALSE) %>%
  mutate(intensity_log2 = log2(fg_quantity)) %>%
  normalise(
    sample = r_file_name,
    intensity_log2 = intensity_log2,
    method = "median"
  ) %>%
  filter(pep_is_proteotypic == TRUE)

data_normalised <- data_normalised %>%
  filter(normalised_intensity_log2 > 5)


# Fit without normalization

fit <- data_normalised %>%
  fit_drc_4p(
    sample = r_file_name,
    grouping = eg_precursor_id,
    response = normalised_intensity_log2,
    dose = r_condition,
    filter = "post",
    n_replicate_completeness = 2, 
    n_condition_completeness = 4, 
    retain_columns = c(pg_protein_accessions)
  )

drc_4p_plot(
  fit,
  grouping = eg_precursor_id,
  dose = r_condition,
  response = normalised_intensity_log2,
  targets = "_VFDVELLKLE_.2",
  unit = "pM",
  x_axis_limits = c(10, 1e+08),
  export = FALSE
)
# make sure to retain columns that you need later but that are not part of the function


# ------------------------------------------------------------------------------
# -- Impute data with Random Forest --------------------------------------------
# ------------------------------------------------------------------------------

library(missForest)

data_wide <- data_normalised %>% 
  tidyr::pivot_wider(
    id_cols=r_file_name,
    names_from = eg_precursor_id, 
    values_from = normalised_intensity_log2)

rownames(data_wide) <- data_wide$r_file_name
file_names <- data_wide$r_file_name

data_wide <- data_wide[, !names(data_wide) %in% c("r_file_name")]
data_wide <- data.frame(lapply(data_wide, function(x) as.numeric(as.character(x))))

data_imputed_rf <- missForest(data_wide, verbose  = TRUE)

data_imputed <- data_imputed_rf$ximp
data_imputed["r_file_name"] <- file_names


data_imputed <- data_imputed %>%
  as.data.frame() %>%                      # Convert to dataframe
  tidyr::pivot_longer(
    cols = -r_file_name,                    # Pivot all columns except r_file_name
    names_to = "eg_precursor_id",           # New column for precursor ID
    values_to = "normalised_intensity_log2" # New column for intensity values
  )

data_normalised <- data_normalised %>%
  mutate(eg_precursor_id = paste0("X", eg_precursor_id))

data_joined <- data_imputed %>%
  full_join(
    data_normalised,
    by = c("r_file_name", "eg_precursor_id")
  )

file_condition_mapping <- data_joined %>%
  filter(!is.na(r_condition)) %>%          # Keep only rows with non-NA r_condition
  select(r_file_name, r_condition) %>%     # Select relevant columns
  distinct()  

data_joined <- data_joined %>%
  left_join(file_condition_mapping, by = "r_file_name") %>%
  mutate(r_condition = coalesce(r_condition.x, r_condition.y)) %>%  # Fill missing r_condition
  select(-r_condition.x, -r_condition.y)  

fit_plot <- drc_4p_plot(
  fit,
  grouping = eg_precursor_id,
  dose = r_condition,
  response = normalised_intensity_log2,
  targets = "_TPVTWLVGPFAPGITEK_.2",
  unit = "pM",
  x_axis_limits = c(10, 1e+08),
  export = FALSE
)

filtered_data <- data_joined %>%
  dplyr::filter(eg_precursor_id == "X_TPVTWLVGPFAPGITEK_.2") %>%
  dplyr::filter(is.na(normalised_intensity_log2.y)) %>%
  dplyr::rename(normalised_intensity_log2 = normalised_intensity_log2.x) %>%
  dplyr::mutate(
    r_condition = as.numeric(r_condition),
    normalised_intensity_log2 = as.numeric(normalised_intensity_log2)
  ) 
  

# Add the filtered red dots to the existing ggplot
plot1 <- fit_plot[[1]]

plot1 + 
  ggplot2::geom_point(
    data = filtered_data,
    ggplot2::aes(x = r_condition, y = normalised_intensity_log2),
    color = "red",   # Red color for the new dots
    size = 3         # Adjust the size as needed
  )





