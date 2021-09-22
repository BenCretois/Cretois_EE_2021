########################################################################
#################### Function to fit the RSF models ####################
########################################################################

# This function fits all RSF models
# - An RSF model made with the telemetry data
# - A naive RSF made with Citizen Science obs & random availability points
# - A corrected RSF made with Citizen Science obs & corrected availability points

# Models use a Weighted logistic regression as suggested by Filthian 2013 & Muff 2019. 
# In their paper, the authors also advise using random effects for all covariates to 
# account for individual variations. This was done for the telemetry data.

# We set the weight of the pseudo-absences to 1000.

fit_rsf_models <- function(data){

  #################################################
  # Fit a RSF for the telemetry data : the truth! #
  #################################################
  inla.setOption(enable.inla.argument.weights=TRUE)
  
  # Get the telemetry obs and the random telemetry
  data_telem <- data %>% 
    filter(type == 'telem_gps' | type == 'random_telem') %>%
    mutate(ID = as.numeric(as.factor(FK_RegNr))) %>% 
    mutate(ID2 = ID, ID3 = ID, ID4 = ID, ID5 = ID, ID6= ID, ID7 = ID, ID8 = ID) %>% 
    mutate(weight = 1000^(1 - value))
  
  ### Model for RSF with telemetry data
  # We fit a model which has both random intercept and random slopes as
  # suggested in Muff et al. 2020
  
  formula_telem1 <- value ~ 1 + d_roads_log + d_urb_log + path_use_log + slope + n_forest_log + n_agr_log + alt_log +
    f(ID, model = 'iid', hyper=list(theta = list(initial=log(1e-6),fixed=T)))+
    f(ID2, d_roads_log, model="iid",hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) + 
    f(ID3, d_urb_log, model="iid",hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) + 
    f(ID4, path_use_log,model="iid",hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) +
    f(ID5, slope,model="iid",hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) + 
    f(ID6, n_forest_log,model="iid",hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) + 
    f(ID7, n_agr_log,model="iid",hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) +
    f(ID8, alt_log,model="iid",hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) 
  
  mod_telem1 <- inla(formula_telem1, 
                           family = 'binomial',
                           data = data_telem, weights = data_telem$weight,
                           verbose = TRUE,
                           control.compute = list(config = TRUE, # Draw samples from posterior
                                                  waic = TRUE, dic = TRUE, cpo = TRUE),
                           control.inla = list(strategy = "gaussian", 
                                               int.strategy = "eb",
                                               reordering = 'amdc'))
  
  
  #########################################
  # Fit a RSF with the "naive" CS dataset #
  #########################################

  # Get the CS observations and the random CS
  data_naive_cs <- data %>% 
    filter(type == 'cs' | type == 'random_cs') %>% 
    mutate(weight =  1000^(1 - value))
  
  # Define the formula for the CS dataset
  formula_cs <- value ~ 1 + d_roads_log + d_urb_log + path_use_log + slope + n_forest_log + n_agr_log + alt_log
  
  # Model
  mod_naive_cs <- inla(formula_cs, 
                             family = 'binomial',
                             data = data_naive_cs, weights = data_naive_cs$weight,
                             verbose = TRUE,
                             control.compute = list(config = TRUE, # Draw samples from posterior
                                                    waic = TRUE, dic = TRUE, cpo = TRUE),
                             control.inla = list(strategy = "gaussian", 
                                                 int.strategy = "eb",
                                                 reordering = 'metis'))
  
  #############################################
  # Fit a RSF with the "corrected" CS dataset #
  #############################################
  
  # Get the CS observations and the random CS
  data_corrected_cs <- data %>% 
    filter(type == 'cs' | type == 'corrected_cs') %>% 
    mutate(weight =  1000^(1 - value))
  
  # Model
  mod_corrected_cs <- inla(formula_cs, 
                                 family = 'binomial',
                                 data = data_corrected_cs, weights = data_corrected_cs$weight,
                                 verbose = TRUE,
                                 control.compute = list(config = TRUE, # Draw samples from posterior
                                                        waic = TRUE, dic = TRUE, cpo = TRUE),
                                 control.inla = list(strategy = "gaussian", 
                                                     int.strategy = "eb",
                                                     reordering = 'metis'))
  
  ########################################################
  # Fit the RSF models with the cross-species correction #
  ########################################################
  
  # Get the CS observations and the corrected availability for species 1
  data_corrected_cross_sp1 <- data %>% 
    filter(type == 'cs' | type == 'corrected_1') %>% 
    mutate(weight =  1000^(1 - value))
  
  # Model
  mod_corrected_cross_sp1 <- inla(formula_cs, 
                           family = 'binomial',
                           data = data_corrected_cross_sp1, weights = data_corrected_cross_sp1$weight,
                           verbose = TRUE,
                           control.compute = list(config = TRUE, # Draw samples from posterior
                                                  waic = TRUE, dic = TRUE, cpo = TRUE),
                           control.inla = list(strategy = "gaussian", 
                                               int.strategy = "eb",
                                               reordering = 'metis'))
  
  # Get the CS observations and the corrected availability for species 2
  data_corrected_cross_sp2 <- data %>% 
    filter(type == 'cs' | type == 'corrected_2') %>% 
    mutate(weight =  1000^(1 - value))
  
  # Model
  mod_corrected_cross_sp2 <- inla(formula_cs, 
                           family = 'binomial',
                           data = data_corrected_cross_sp2, weights = data_corrected_cross_sp2$weight,
                           verbose = TRUE,
                           control.compute = list(config = TRUE, # Draw samples from posterior
                                                  waic = TRUE, dic = TRUE, cpo = TRUE),
                           control.inla = list(strategy = "gaussian", 
                                               int.strategy = "eb",
                                               reordering = 'metis'))
  
  ##################
  # Return objects #
  ##################
  plot_coef1 <- Efxplot(list(mod_telem1, mod_naive_cs, mod_corrected_cs),
                       ModelNames = c('Telemetry RS', 'Naive CS', 'Corrected CS')) +
    coord_flip(ylim = c(-0.5,0.5))
  
  plot_coef2 <- Efxplot(list(mod_telem1, mod_naive_cs, mod_corrected_cs, mod_corrected_cross_sp1, mod_corrected_cross_sp2),
                        ModelNames = c('Telemetry RS', 'Naive CS', 'Corrected CS', 'Corrected cross sp 1', 'Corrected cross sp 2')) +
    coord_flip(ylim = c(-0.5,0.5))
  
  return(list(plot_coef1, plot_coef2, mod_telem1, mod_naive_cs, 
              mod_corrected_cs, mod_corrected_cross_sp1, mod_corrected_cross_sp2))
}
