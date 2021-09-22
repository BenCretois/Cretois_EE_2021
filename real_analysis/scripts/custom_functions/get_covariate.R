##########################################################################
#################### Function to get covariate values ####################
##########################################################################

# This script aims to extract the covariate values for each observation. 
# The n argument is useful for setting where, in the dataset, to begin the
# extraction. For instance

get_covariate <- function(species, data, list_raster, n){

  data <- data %>% 
    mutate(d_roads = NA, d_urb = NA, altitude = NA, pop = NA, n_agr = NA, n_forest = NA, path_use = NA, slope = NA)
    
  if(species == 'wreindeer'){

    for(i in 1:length(list_raster)){
      r <- read_stars(paste0("real_analysis/covariates_reindeer/", list_raster[[i]]))
      data[ , i+n] <- raster_extract(r, data)
      rm(r)
    }
  }
    
  else{

      for(i in 1:length(list_raster)){
        r <- read_stars(paste0("real_analysis/covariates_moose_roe/", list_raster[[i]]))
        data[ , i+n] <- raster_extract(r, data)
        rm(r)
      }
    }
    
    # Somehow distance variables are lists inside the df -> need to vectorize them
    data$d_roads = as.vector(unlist(data$d_roads))
    data$d_urb = as.vector(unlist(data$d_urb))
    data$slope = as.vector(unlist(data$slope))
    data$n_forest = as.vector(unlist(data$n_forest))
    data$path_use = as.vector(unlist(data$path_use))
    data$pop = as.vector(unlist(data$pop))
    data$altitude = as.vector(unlist(data$altitude))
    
    # Replace the NAs by 0 for pop, n_agr and n_forest
    data <- data %>% replace_na(list(pop = 0, n_agr = 0, n_forest = 0))
    
    # Log the variables
    data <- data %>% mutate(
      d_roads_log = log1p(d_roads),
      d_urb_log = log1p(d_urb),
      n_agr_log = log1p(n_agr),
      n_forest_log = log1p(n_forest),
      path_use_log = log1p(path_use),
      pop_log = log1p(pop),
      alt_log = log1p(altitude))
  
  # Return the full dataset
  return(data)
}

