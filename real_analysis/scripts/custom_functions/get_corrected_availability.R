##################################################################################
# Function computing the bias model and availability points for moose & roe deer #
##################################################################################

# We want to estimate how different are CS observations from the telemetry data. 
# We fit a binomial "observer model" where the 0s are the telemetry points while 
# the 1s are the Citizen Science observations.

# After fitting the observer model we used randomly sampled points to compute the 
# probability of having a CS observations at a specific location, this represents 
# the corrected availability for Citizen Scientists. We sample a number of corrected 
# availability for fitting the "corrected CS model" later in the analysis.

get_corrected_availability <- function(species, data, formula_bias, ratio_availability, list_raster){

  #####################
  # Dataset wrangling #
  #####################
  
  data <- data %>% 
    filter(type == 'telem_gps' | type == 'cs') %>% 
    mutate(type_recode = ifelse(type == 'cs', 1, 0)) # Here we code CS obs as 1 and telemetry as 0
  
  ###################################
  # Get data for spatial prediction #
  ###################################
  
  # Each species has multiple study sites, we need to take each study site and not
  # the entire area -> draw multiple Convex Polygons
  
  # Wrangle to get the minimum convex polygon for each animal
    for_mcp <- data %>% 
      mutate(id = as.factor(FK_RegNr)) %>% 
      dplyr::select(id) %>% 
      as_Spatial()
    
    mcp <- mcp(for_mcp, percent = 97)

    # Merge the polygons
    mcp_sf <- st_as_sf(mcp)
    mcp_union <- st_union(mcp_sf)
    mcp <- st_buffer(mcp_union, dist = 10000) %>% 
      st_as_sf(crs = 25833)
  
    # Sample random points in this area
  s <- st_sample(mcp, 100000) %>% st_as_sf(crs = 25833)
  
  coords_s <- st_coordinates(s) %>% as.tibble() 
  
  s <- s %>% st_drop_geometry() %>% 
    bind_cols(., coords_s) %>% 
    st_as_sf(coords = c('X','Y'), crs = 25833) 
  
  ############################
  # Extract covariate values #
  ############################
  
  if(species == 'wreindeer'){
    # Extract the covariates for the prediction points
    s <- get_covariate('wreindeer', s, list_raster, 1)
  }
  
  else{
    s <- get_covariate('other', s, list_raster, 1)
  }
  
  # Need to drop NA for path use as some of the points are outside the raster dimensions
  s <- s %>% 
    drop_na()
  
  ################################
  # Prepare dataset for analysis #
  ################################
  
  # log the variables -> all very skewed
  data_prep <- data %>% mutate(
    weight = 1000^(1 - type_recode)) %>% 
    st_drop_geometry()
  
  if(species == 'wreindeer'){
    data_prep <- sample_n(data_prep, 40000)
  }
  else{
    data_prep <- data_prep
  }
  
  # Log the predictors for the prediction points too
  s_prep <- s %>%
    dplyr::select(path_use_log, d_roads_log, d_urb_log, n_forest_log, n_agr_log, alt_log, slope, pop_log) %>% 
    st_drop_geometry()
  
  # Take the right part of the formula bias
  f <- splitFormula(formula_bias, sep = '~')[[1]]
  
  # Create the model matrix
  Xmat <- model.matrix(f, data = s_prep) %>% as.data.frame()
  
  ######################
  # Fit the bias model #
  ######################
  
  # Weighted logistic regression as suggested by Filthian 2013 & Muff 2019
  inla.setOption(enable.inla.argument.weights=TRUE)
  
  # Prior on the fixed effects
  prior.fixed <- list(mean.intercept = 0, prec.intercept = 1,
                      mean = 0, prec = 1)
  
  formula <- formula_bias

  bias_model <- inla(formula, 
                     family = 'binomial',
                     data = data_prep, 
                     weights = data_prep$weight,
                     control.inla = list(strategy = "gaussian", 
                                         int.strategy = "eb", 
                                         reordering = 'metis'),
                     verbose = TRUE)

  ###################
  # Check the model #
  ###################
  
  # Prepare the dataset for checks
  #data_prep <- data_prep %>% mutate(fitted = bias_model$summary.fitted.values$mean,
   #                                 observed = data$type_recode) %>% 
  #  mutate(resid = observed - fitted, logit = log(fitted / (1 - fitted))) %>% 
  #  mutate(pearson_resid = (resid - mean(resid))^2 / resid, ID = 1:nrow(.))
  
  #check_plot <- data_prep %>% 
  #  dplyr::select(logit, path_use_log, d_roads_log, d_urb_log, pop_log, d_roads_log, n_forest_log) %>% 
  #  gather(key = 'predictors', value = 'predictors_value', -logit)
  
  # Check if the log-odds are linearly related to the predictors
  #check1 <- ggplot(check_plot, aes(y = logit, x = predictors_value)) +
  #  geom_point(size = 0.5, alpha = 0.5) +
  #  geom_smooth() + 
  #  theme_bw() + 
  #  facet_wrap(~predictors, scales = "free_y") +
  #  xlab('Probability of a CS observation')

  # Plot Pearson's residuals values
  #check2 <- ggplot(data_prep, aes(ID, pearson_resid)) + 
  #  geom_point(aes(color = observed), alpha = .5) +
  #  theme_bw()
  
  ##############################################################################
  # Predict the probability of having a CS observation at the sampled location #
  ##############################################################################
  
  # Make the predictions
  X_mat <- as.matrix(Xmat)
  coef_mat <- as.matrix(bias_model$summary.fixed$mean)
  s$pred <- inv.logit(X_mat %*% coef_mat) 
  
  ########################################################
  # Sample random points with regards to the probability #
  ########################################################
  
  cs_availability <- data.frame(type = NA, value = NA)
  n_sample <- data %>% filter(type == 'cs') %>% nrow(.)
  n_sample <- ratio_availability * n_sample
  
  cs_availability <- s %>% 
    sample_n(., size = n_sample, weight = s$pred, replace = TRUE) %>% 
    mutate(type = 'corrected_cs', value = 0) %>% 
    dplyr::select(type, value) %>% 
    st_as_sf(crs = 25833)
  
  ##################
  # Return objects #
  ##################
  
  # The function return:
  # - the the corrected availibility points
  # - the bias model
  return(list(cs_availability, bias_model))
}

