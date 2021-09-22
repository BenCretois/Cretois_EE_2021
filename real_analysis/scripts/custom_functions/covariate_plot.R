################################################################################
# Script to plot the covariate distribution for CS, telemetry and availability #
################################################################################

# The function returns a boxplot of the distribution of the covariates for each 
# "type" of observations: Telemetry, Citizen Science and availability

covariate_plot <- function(data){
  
  # Take only telem, CS and availability
  data <- data %>% 
    filter(type == 'telem_gps' | type == 'cs' | type == 'random_cs')
  
  # Make the plot
  data_descriptive <- data %>% st_drop_geometry() %>%  
    dplyr::select(type, path_use_log, d_roads, d_urb, n_forest, n_agr, slope, altitude) %>% 
    gather(key = 'predictors', value = 'predictors_value', -type)
  
  # plot boxplots of all variables for all "types"
  cov_plots_tot <- ggplot(data_descriptive, aes(x = type, y = predictors_value, fill = type)) +
    geom_boxplot() +
    geom_smooth() + 
    theme_bw() + 
    facet_wrap(~predictors, scales = "free_y") +
    scale_fill_viridis(option = 'D', discrete = TRUE, alpha = .7)
  
  # plot for only distance to roads, urban center and path use intensity
  data_descriptive2 <- data %>% st_drop_geometry() %>% filter(type != 'random_telem') %>% 
    dplyr::select(type, path_use_log, d_roads, d_urb) %>% 
    gather(key = 'predictors', value = 'predictors_value', -type) 
  
  new_labels <- c('d_roads'='dist. to roads', 'd_urb'='dist. to urb. settlements', 'path_use_log'='log path use intensity')
  
  ggplot(data_descriptive2, aes(x = predictors_value, y = type, fill = type)) +
    geom_boxplot() +
    #geom_density_ridges() +
    geom_smooth() + 
    theme_bw() + 
    facet_wrap(~predictors, scales = "free", labeller = labeller(predictors = new_labels)) +
    ylab(' ') + xlab('Predictor value') + scale_fill_manual(name = "Dataset", 
                                                            labels = c("Citizen Science", "Availability", "Telemetry"),
                                                            values = c("#f7cb44ff", "#de7065ff", "#7e4e90ff")) + 
    scale_y_discrete(labels = c("cs" = "Citizen science", "random_cs" = "availability", "telem_gps" = "Telemetry"))
  
  
  # Object to return
  return(list(cov_plots_tot, cov_plots_paper))
}






