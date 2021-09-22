#############################################
# Script simulation for Citizen Science RSF #
#############################################

# We want to evaluate whether sampling availability
# points with regard to a bias model can improve 
# the estimation of parameters meaningful for species'
# ecology.

# Load libraries
library(NLMR)
library(raster)
library(sf)
library(stars)
library(geobgu) # download from github
library(tidyverse)
library(INLA)
library(ggregplot) # download from github
library(nlme)

##########################
# Small helper functions #
##########################
inv.logit <- function(x) exp(x) / (1 + exp(x))
my_rmse <- function(truth, predicted) sqrt(mean((truth - predicted)^2))

#################################################
# Function to simulate environmental covariates #
#################################################


simul_env <- function(dim){
  # Distance to urban centers
  d_urb <- nlm_distancegradient(ncol = dim[1], nrow = dim[2], resolution = 1, 
                                origin = c(3,5,6,7)) %>% st_as_stars()
  # Presence / absence of forests
  forests <- nlm_randomcluster(ncol = dim[1], nrow = dim[2],p = 0.4, 
                               ai = c(0.40, 0.60)) %>% st_as_stars()
  # A gradient correlated with roads ~ -80
  other_grad <- nlm_edgegradient(dim[1], nrow = dim[2], direction = 160)
  # Altitude gradient
  altitude <- nlm_gaussianfield(dim[1], nrow = dim[2], mean = 1) %>% 
    st_as_stars()
  # Distance to roads
  inv_raster <- function(x) 1 - x
  d_roads <- nlm_edgegradient(dim[1], nrow = dim[2], direction = 180) %>% 
    inv_raster(.) %>% st_as_stars()
  # Presence / absence of a nice viewpoint
  nice_viewpoints <- nlm_randomcluster(ncol = dim[1], nrow = dim[2], p = 0.4, 
                                       ai = c(0.90, 0.10)) %>% st_as_stars()
  # list the covariates
  list_raster <- list(d_urb, forests, other_grad, altitude, d_roads, 
                      nice_viewpoints)
  return(list_raster)
}

#########################################################
############### Simulate the data #######################
#########################################################

simul_data <- function(dim,
                       list_raster,
                        beta0_sp,
                        beta0_cs,
                        beta_forests,
                        beta_other_gradient,
                        beta_altitude,
                        beta_roads_sp,
                        beta_roads_cs,
                        beta_urb,
                        beta_nicev,
                        plotdat = TRUE,
                        free_ditribution = TRUE){
# Set the dimensions
dim = dim

############################
# Set the parameter values #
############################

# For the species
beta0_sp = beta0_sp; beta_forests = beta_forests; beta_other_gradient; beta_altitude = beta_altitude; beta_roads_sp = beta_roads_sp

# For the citizen science
beta0_cs = beta0_cs; beta_roads_cs = beta_roads_cs; beta_urb = beta_urb; beta_nicev = beta_nicev 

#######################################################################
# Sample presence / absence of Citizen Science observations & species #
#######################################################################

# Create a grid and take the centroids
bbox <- st_bbox(list_raster[[1]]) %>% st_as_sfc()
grid <- st_make_grid(bbox, n = dim, what = 'polygons') %>% st_as_sf()
centr <- st_centroid(grid) %>% st_as_sf() %>% 
  mutate(d_urb = NA, forests = NA, other_gradient = NA, altitude = NA, 
         d_roads = NA, nice_viewpoints = NA)

# Extract the covariates for each of the centroids
for(i in 1:length(list_raster)){
  centr[ , i+1] <- raster_extract(list_raster[[i]], centr)
}

# Create the presence / absence dataset 

centr <- centr %>% 
  mutate(proba_species = inv.logit(beta0_sp + beta_forests * forests + 
                                     beta_other_gradient * other_gradient + 
                                     beta_altitude * altitude + 
                                     beta_roads_sp * d_roads)) %>% 
  mutate(species_occ = rbinom(nrow(.), 1, proba_species)) %>% 
  # If free distribution = TRUE then probability of CS obs = P(CS present) 
  # * P(Species use habitat)
  mutate(proba_cs = inv.logit(beta0_cs + beta_roads_cs * d_roads + 
  beta_urb * d_urb + beta_nicev * nice_viewpoints) * proba_species) %>% 
  mutate(cs_occ = rbinom(nrow(.), 1, proba_cs)) 

##################
# Return objects #
##################

if(plotdat == TRUE){
  for_plot <- centr %>% dplyr::select(species_occ, 
                                      cs_occ)
  plot(for_plot)
}
return(centr)
}

#########################################################################################
######### Functions which fit all the models (truth, naive and corrected RSF) ###########
##########################################################################################

fit_models <- function(data, env,
                       formula_bias, size_buffer, return_coefs = FALSE){

  # Set the weights in INLA
  inla.setOption(enable.inla.argument.weights=TRUE)
  
  ##############################
  # Make the different dataset #
  ##############################
  species_pres <- data %>% 
    filter(species_occ == 1) %>% 
    mutate(bias = 0)
  
  cs_pres <- data %>% 
    filter(cs_occ == 1) %>% 
    mutate(bias = 1)
  
  ### For the species model 
  species_full <- data %>% 
    filter(species_occ == 0) %>% 
    sample_n(., size = 2*nrow(species_pres)) %>% 
    bind_rows(., species_pres) %>% 
    mutate(weight = 1000^(1 - species_occ)) 
  
  ### For the naive RSF model
  cs_naive <- data %>% 
    filter(cs_occ == 0) %>% 
    sample_n(., size = 5*nrow(cs_pres)) %>% 
    bind_rows(., cs_pres) %>% 
    mutate(weight = 1000^(1 - cs_occ)) 
  
  ### For the bias model
  df_biases <- bind_rows(species_pres, cs_pres) %>% 
    mutate(weight = 1000^(1 - bias)) 
  
  ##################
  # Fit the models #
  ##################
  
  # True RSF
  mod_sp <- inla(species_occ ~ forests + altitude + d_roads + other_gradient, 
                 data = species_full, 
                 family = 'binomial', weights = species_full$weight,
                 control.inla = list(strategy = "gaussian", 
                                     int.strategy = "eb",
                                     reordering = 'amdc'))
  
  # Bias model
  mod_bias <- inla(formula_bias, 
                   data = df_biases, 
                   family = 'binomial', weights = df_biases$weight,
                   control.inla = list(strategy = "gaussian", 
                                       int.strategy = "eb",
                                       reordering = 'amdc'))
  
  # Naive RSF with uncorrected CS observations
  rsf_naive <- inla(cs_occ ~ forests + altitude + d_roads + other_gradient, 
                    data = cs_naive, 
                    family = 'binomial', weights = cs_naive$weight,
                    control.inla = list(strategy = "gaussian", 
                                        int.strategy = "eb",
                                        reordering = 'amdc'))
  
  #################
  # Corrected RSF #
  #################
  
  # Randomly sample points across the landscapes to compute the probability they host a CS observation
  bbox <- st_bbox(data) %>% st_as_sfc()
  random_points <- st_sample(bbox, 100000) %>% st_sf() %>% 
    mutate(d_urb = NA, forests = NA, other_gradient = NA, altitude = NA, 
           d_roads = NA, nice_viewpoints = NA)
  
  # Get the covariate values at each of these points
  for(i in 1:length(env)){
    random_points[ , i+1] <- raster_extract(env[[i]], random_points)
  }

  # We sample availability point from the bias model
  pred_bias_f <- splitFormula(formula_bias, sep = '~')[[1]]
  Xmat <- model.matrix(pred_bias_f, data = random_points)
  
  X_mat <- as.matrix(Xmat)
  coef_mat <- as.matrix(mod_bias$summary.fixed$mean)
  
  random_points <- random_points %>%  mutate(pred_bias = inv.logit(X_mat %*% coef_mat))
  random_points <- sample_n(random_points, 2000, weight = random_points$pred_bias)

  # We create the dataset
  cs_corrected <- random_points %>% 
    mutate(cs_occ = 0) %>% 
    bind_rows(., cs_pres) %>% 
    mutate(weight = 1000^(1 - cs_occ))
  
  rsf_corrected_df <- bind_rows(cs_pres, cs_corrected) %>% 
    mutate(weight = ifelse(cs_occ == 1, 1, 1000))
  
  # Fit the model
  rsf_corrected <- inla(cs_occ ~ forests + altitude + d_roads + other_gradient, 
                        data = rsf_corrected_df, 
                        family = 'binomial', weights = rsf_corrected_df$weight,
                        control.inla = list(strategy = "gaussian", 
                                            int.strategy = "eb",
                                            reordering = 'amdc'))
  
  ##########################
  # Objects to be returned #
  ##########################
  Xmat_pred <- model.matrix(~ forests + altitude + d_roads + other_gradient, 
                            data = data)
  pred_truth <- inv.logit(Xmat_pred %*% mod_sp$summary.fixed$mean)
  pred_naive <- inv.logit(Xmat_pred %*% rsf_naive$summary.fixed$mean)
  pred_corrected <- inv.logit(Xmat_pred %*% rsf_corrected$summary.fixed$mean)
  n <- nrow(data)
  
  # Calculate the RMSE
  rmse_truth_naive <- my_rmse(pred_truth, pred_naive)
  rmse_truth_corrected <- my_rmse(pred_truth, pred_corrected)
  vec_rmse <- c(rmse_truth_naive, rmse_truth_corrected)
  
  # Compare the coefficients of the 3 models
  p_coef <- Efxplot(list(mod_sp, rsf_naive, rsf_corrected), 
                    ModelNames = c('Truth', 'RSF naive', 'RSF corrected')) 
  p_coef
  # Check the number of availability points I have
  n_ava_naive <- cs_naive %>% filter(cs_occ == 0) %>% st_drop_geometry %>% 
    summarise(n = n())
  n_ava_corrected <- cs_corrected %>% filter(cs_occ == 0) %>% 
    st_drop_geometry %>% summarise(n = n())
  vec_ava <- c(n_ava_naive, n_ava_corrected)
  
  if(return_coefs == TRUE){
  # calculate direction of coef changes - is it an improvement?
  improvement <- rep(NA, length(as.matrix(mod_sp$summary.fixed$mean)))
  
  test1 <- as.matrix(mod_sp$summary.fixed$mean) - as.matrix(rsf_naive$summary.fixed$mean)
  test2 <- as.matrix(mod_sp$summary.fixed$mean) - as.matrix(rsf_corrected$summary.fixed$mean)
  
  improvement <- abs(test1)>abs(test2)
  
  # test3 - does it improve the sign?
  test3 <- sign(as.matrix(mod_sp$summary.fixed$mean))==sign(as.matrix(rsf_naive$summary.fixed$mean))
  if(length(which(test3 == FALSE))>0){
   improvement[(which(test3 == FALSE))] <- sign(as.matrix(mod_sp$summary.fixed$mean[(which(test3 == FALSE))]))==sign(as.matrix(rsf_corrected$summary.fixed$mean[(which(test3 == FALSE))]))
  }}

# Return
if(return_coefs == FALSE){return(list(p_coef, vec_rmse, vec_ava))}
if(return_coefs == TRUE){return(list(p_coef, vec_rmse, vec_ava, improvement))}
}

##################################
# Function to bootstrap the RMSE #
##################################

boot_rmse <- function(data, env,
                      formula, boot, return_coefs = FALSE){
  
list_naive <- list()
list_corrected <- list()
list_improved <- list()

for(i in 1:boot){
  if(return_coefs==TRUE){mod <- fit_models(data, env, formula, return_coefs=TRUE)}else{
  mod <- fit_models(data, env, formula, size_buffer = 10)}
  v_rmse <- mod[[2]]
  list_naive[[i]] = v_rmse[1]
  list_corrected[[i]] = v_rmse[2]
  if(return_coefs == TRUE){list_improved[[i]] = mod[[4]]}
}

df_rmse <- tibble(Naive = unlist(list_naive),
                  Corrected = unlist(list_corrected))

df_difference <- tibble(diff = df_rmse$Naive - df_rmse$Corrected) 

df_mean_sd <- tibble(mean = mean(df_difference$diff),
                     sd = sd(df_difference$diff))

if(return_coefs == TRUE){
improvement <- as.data.frame(list_improved)

improvement_perc <- (rowSums(improvement)/ncol(improvement))*100}

plot <- df_rmse %>%
  gather(., key =  'Model', value = 'rmse') %>% 
  ggplot(., aes(x = rmse, fill = Model)) +
  geom_histogram(alpha = .5) +
  scale_fill_manual(values = c('darkgreen', 'blue'))

# Objects to return
if(return_coefs == FALSE){return(list(df_mean_sd, plot))}
if(return_coefs == TRUE){return(list(df_mean_sd, improvement_perc))}
}

#########################################################################################
######### Function to run a full analysis ###########
##########################################################################################


# function to run a full analysis to allow bootstrap of the data generation
# process as well as running the model

# INPUTS: dimensions of domain, formula for the bias model, number of bootstrap for
# selection of background points

run_analysis <- function(dimensions = c(100,100),
                         formula_bias,
                       boot){

# Set the dimensions for the simulation
dim <- dimensions

# Simulate environmental variables
env <- simul_env(dim)

# Simulate data
data <- simul_data(dim,
                   env,
                   beta0_sp = -7, # intercept for SP
                   beta0_cs = 2, # intercept for CS
                   beta_forests = 2.5, # effect of forest cover
                   beta_other_gradient = 3, # effect of the gradient correlated with roads
                   beta_altitude = -2, # effect of altitude
                   beta_roads_sp = 4.5, # effect of distance to road on SP
                   beta_roads_cs = -6, # effect of distance to road on CS
                   beta_urb = -3, # effect of distance to urban centre
                   beta_nicev = 1, # effect of nice viewpoints
                   plotdat = TRUE)

# Run the model
model <- fit_models(data, env, formula_bias, return_coefs = TRUE)

# Run bootstrap of the background point selection
if(boot > 1){bootstrapped_model <- boot_rmse(data, env, 
                                             formula_bias, 
                                             boot, 
                                             return_coefs = TRUE)}

# OUTPUTS: 
# save out the model and the bootstrap as a list

if(boot > 1){output <- list(model, bootstrapped_model)}else{output <- list(model)}

return(output)

}

