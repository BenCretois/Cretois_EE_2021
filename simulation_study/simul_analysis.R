#########################################
# Script run simulations for appendix #
#########################################

# aim of this code:

# generate results bootstrapping the data generation (n = 50)
# generate results bootstrapping the model (n = 50)
# generate results double bootstrapping (n = 30 each)

# present all results in terms of covariate improvement 
# AND RMSE

#########################################

source('simulation_study/simul_functions.R')

library(patchwork)
library(purrr)
library(ggplot2)

# Set the dimensions for the simulation
dimensions = c(100,100)

# Set bias formula

formula_bias <- bias ~ d_roads + d_urb + nice_viewpoints

# run an analysis for the bootstrapping the data generation process
# using 'replicate'

bootstrap_data_gen_50 <- replicate(100, run_analysis(dimensions = dimensions,
                                                   formula_bias = formula_bias,
                                                   boot = 1))

save(bootstrap_data_gen_50, file = "bootstrap_data_gen_50.RData")

# run an analysis for the bootstrapping the model

bootstrap_model_50 <- replicate(1, run_analysis(dimensions = dimensions,
                                                   formula_bias = formula_bias,
                                                   boot = 3))

#save(bootstrap_model_50, file = "bootstrap_model_50.RData")

#bootstrap_both_30 <- replicate(30, run_analysis(dimensions = dimensions,
                                                   formula_bias = formula_bias,
                                                   boot = 30))

#save(bootstrap_both_30, file = "bootstrap_both_30.RData")

# Last one is too slow for now

#########################################
# Plots
#########################################

load("bootstrap_data_gen_50.RData")
load("bootstrap_model_50.RData")

#### Proportion of times covariate is correctly changed ####

# first collapse to the data I want

bootstrap_data <- map_df(bootstrap_data_gen_50, ~{
  row <- data.frame(forests = .x[[4]][1],
                    altitude = .x[[4]][2],
                    d_roads = .x[[4]][3],
                    other_gradient = .x[[4]][4])
  return(row)
  })

# then calculate proportions

coef_results_data <- data.frame(proportion = colSums(bootstrap_data)/100,
                          covariate = colnames(bootstrap_data))

# and plot

coef_plot_data <- ggplot(aes(x=covariate, y=proportion), 
                    data = coef_results_data) +
  geom_col(position = "dodge2") +
  theme_minimal()

coef_plot_data

ggsave(last_plot(), filename = "bootstrap_of_data.png")

#### RMSE and proportion of times covariate is correctly changed ####

# first collapse to the data I want

coef_results_model <- bootstrap_model_50 %>% 
  lapply(`[[`, 2) 

coef_results_model <- coef_results_model[[2]] %>%
  setNames(c("Int", "Forests", "Altitude", 
             "D_roads", "Other_gradient")) %>%
  as.data.frame()

coef_results_model <- data.frame(covariate = rownames(coef_results_model),
                                proportion = coef_results_model[,1]) %>%
  filter(covariate != "Int")

# and plot

coef_plot_model <- ggplot(aes(x=covariate, y=proportion), 
                    data = coef_results_model) +
  geom_col(position = "dodge2") +
  theme_minimal()

coef_plot_model

ggsave(last_plot(), filename = "bootstrap_of_model.png")




