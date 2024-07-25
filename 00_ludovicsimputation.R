# ludo imputation

##### imputation functions with fixed set.seeds
impute_MNAR <- function(n_missing, min, sd, rnk){
  rnorm(n_missing, mean = min - log2(3), sd = sd)[rnk]
}

impute_MAR <- function(n_missing, mean, sd, rnk){
  rnorm(n_missing, mean = mean, sd = sd)[rnk]
}


impute_ludo <- function(data, r_condition, r_condrep, precursor, intensity){
  data_imputed <- data %>%
    dplyr::select({{r_condrep}}, {{r_condition}}, {{precursor}}, {{intensity}})%>%
    tidyr::complete({{precursor}}) %>%
    dplyr::arrange(!is.na({{intensity}}), {{r_condrep}}) %>% 
    group_by({{precursor}}, {{r_condition}})%>%
    dplyr::mutate(mean = mean({{intensity}}, na.rm = TRUE)) %>%
    dplyr::mutate(sd = sd({{intensity}}, na.rm = TRUE))%>%
    dplyr::mutate(missing = sum(is.na({{intensity}})))%>%
    dplyr::mutate(repl = dplyr::n())%>%
    dplyr::mutate(rnk = 1:dplyr::n())%>%
    group_by({{precursor}})%>%
    dplyr::mutate(sd = mean(sd, na.rm = TRUE))%>%
    dplyr::mutate(min = min({{intensity}}, na.rm = TRUE))%>%
    rowwise() %>%
    dplyr::mutate(NA_type = ifelse(is.na({{intensity}}),
                                   ifelse(repl - missing == 0, "MNAR", "MV"), "complete")) %>%
    dplyr::mutate(normalised_intensity_imputed_log2 = ifelse(is.na({{intensity}}),
                                                             ifelse(repl - missing == 0, impute_MNAR(missing, min, sd, rnk), impute_MAR(missing, mean, sd, rnk)), {{intensity}})) %>%
    dplyr::mutate(imputed = is.na({{intensity}}))%>%
    dplyr::select({{r_condrep}}, {{precursor}}, normalised_intensity_imputed_log2, NA_type, imputed)%>%
    ungroup()
  
  data %>%
    left_join(data_imputed, by = c(as_label(enquo(r_condrep)), as_label(enquo(precursor))))
}
