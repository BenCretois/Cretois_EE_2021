#############################################
# Function to run a full analysis #
#############################################

# function to run a full analysis to allow bootstrap of the data generation
# process as well as running the model

# INPUTS: dimensions of domain, formula for the bias model, number of bootstrap for
# selection of background points

run_analysis <- function(dim = c(100,100),
                         formula_bias,
                         boot){
  
# source scripts and load required packages
source('simul_functions.R')
  
# Set the dimensions for the simulation
dim = dim
  
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
model <- fit_models(data, formula_bias, return_coefs = TRUE)

# Run bootstrap of the background point selection
if(boot > 1){bootstrapped_model <- boot_rmse(data, 
                                formula_bias, 
                                boot, 
                                return_coefs = TRUE)}

# OUTPUTS: 
# save out the model and the bootstrap as a list

if(boot > 1){output <- list(model, bootstrapped_model)}else{output <- list(model)}

return(output)

}