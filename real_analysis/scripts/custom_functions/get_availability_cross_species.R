################################################################################
# Script to get corrected availability using the bias model of another species #
################################################################################

get_availability_cross_species <- function(species, data, species_bias_model, bias_model, ratio_availability, list_raster, formula_bias, type){
  
  # Get the MCP for each animal
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
  
  if(species_bias_model == 'wreindeer'){
    s_prep <- s %>%
      dplyr::select(path_use_log, d_roads_log, d_urb_log, n_forest_log, n_agr_log, alt_log, slope) %>% 
      st_drop_geometry()
  }
  else{
    s_prep <- s %>%
      dplyr::select(path_use_log, d_roads_log, d_urb_log, n_forest_log, pop_log, n_agr_log, alt_log, slope) %>% 
      st_drop_geometry()
  }
  
  ##############################################################################
  # Predict the probability of having a CS observation at the sampled location #
  ##############################################################################
  
  # Take the right part of the formula bias
  f <- splitFormula(formula_bias, sep = '~')[[1]]
  
  # Create the model matrix
  Xmat <- model.matrix(f, data = s_prep) %>% as.data.frame()
  
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
    mutate(type = type, value = 0) %>% 
    dplyr::select(type, value) %>% 
    st_as_sf(crs = 25833)
  
  ##################
  # Return objects #
  ##################
  
  return(cs_availability)  
}
