###################################################################
##################### Make the Moose dataset ######################
###################################################################

make_moose_dataset <- function(buffer_for_cs, buffer_for_telem, thinning, ratio_random_cs){
  
  ######################################
  # Wrangle the telemetry observations #
  ######################################
  
  # Open the dataset and make the necessary transformations:
  # Change CRS
  # Make a "data" column readable by the amt package
  # Add a season index -> we will just need summer in our analysis
  # Open telemetry
  moose_telem <- read_sf("real_analysis/dataset/moose/gps_data_moose.shp") %>% 
    st_transform(crs = 25833) %>% 
    mutate(t = as.POSIXct(acquisitio),
           id = as.factor(animal_id),
           FK_RegNr = animal_id,
           value = 1,
           type = 'telem_gps',
           X = st_coordinates(.[,1]),
           Y = st_coordinates(.[,2])) %>% 
    mutate(date = t) %>% 
    cutData(., type = 'season') %>% 
    dplyr::filter(season == "summer (JJA)") %>% 
    mutate(H = as.numeric(format(t, "%H"))) %>% 
    filter(H > 8 & H < 22) %>% 
    dplyr::select(id, FK_RegNr, t, value, type, season, H, X, Y)


# Define the study area for cropping outliers (points in the ocean)
  bbox <- read_sf("real_analysis/dataset/study_area/study_area.shp") %>% 
    st_transform(crs = 25833)
  moose_crop <- st_intersection(moose_telem, bbox)
  
  # Summarise number of observation per GPS collar // looks good
  id_moose <- moose_crop %>% 
    st_drop_geometry() %>% 
    group_by(FK_RegNr) %>% 
    summarise(n = n()) %>% 
    filter(n > 100)
  
  # Take only the moose with more than 100 relocations
  moose_crop <- moose_crop %>% filter(FK_RegNr %in% id_moose$FK_RegNr)
  

  ###################################################################
  # Thin the telemetry data to avoid spatial / temporal correlation #
  ###################################################################
  vector_id <- unique(moose_crop$FK_RegNr)
  moose_crop_thinned <- list()
  
  for(i in 1:length(vector_id)){
    dat_track <- moose_crop %>% dplyr::filter(FK_RegNr == paste(vector_id[i])) %>% 
    amt::make_track(.x = X, .y = Y, .t = t, key = FK_RegNr)
    dat_thinned <- amt::track_resample(dat_track, rate = thinning, tolerance = hours(1))
    moose_crop_thinned[[i]] <- dat_thinned
  }
  
  moose_crop_thinned <- do.call(rbind, moose_crop_thinned) %>%  
    dplyr::select(FK_RegNr = key, geometry) %>% 
    mutate(value = 1, type = 'telem_gps')
  
  # ---> Check if I have enough telemetry to make the MCP
  
  #############################################################
  # Wrangle to get the minimum convex polygon for each animal #
  #############################################################
  
  moose_sp <- moose_crop_thinned %>% 
    dplyr::select(FK_RegNr, geometry) %>% 
    as_Spatial(.)
  
  mcp <- mcp(moose_sp, percent = 99.5)

    # Merge the polygons and make a buffer around
  # Buffer helps to sample more CS observations,
  # otherwise there is not enough observations
  mcp_sf <- st_as_sf(mcp)
  mcp_union <- st_union(mcp_sf)
  mcp_buffer <- st_buffer(mcp_union, dist = buffer_for_cs) %>% # 50000
    st_as_sf(crs = 25833)
  
 # Crop the observations out of the polygon
  data_telem <- st_intersection(moose_crop_thinned, mcp_buffer)
  
  ######################################
  # Get random points telemetry points #
  ######################################
  
  # For each presence point we sample 2 absence points
  # Enough for the Infinitively Weighted Logistic Regression
  # (Muff et al. 2020)
  
  # Compute the number of points that I need for each animal
  n_points_id <- data_telem %>% 
    st_drop_geometry() %>%
    mutate(id = as.factor(FK_RegNr)) %>% 
    group_by(id) %>% 
    summarise(n = n()) %>% 
    mutate(to_sample = 2*n)
  
  # Join the dataset to get the number of points to sample for each MCP
  mcp_sf <- mcp_sf %>% full_join(., n_points_id, by = 'id') %>% drop_na()
  
  # Make a list of mcps
  randomp <- list()
  
  # Run the sampling algorithm
  for(i in 1:nrow(mcp_sf)){
    mcp <- mcp_sf[i,]
    mcpb <- st_buffer(mcp, buffer_for_telem) # 30000
    p <- mcpb %>% st_sample(., mcpb$to_sample) %>% st_as_sf()
    coords_p <- st_coordinates(p)
    p <- p %>% st_drop_geometry() %>% 
      mutate(X = coords_p[,1],
             Y = coords_p[,2],
             id = mcpb$id,
             value = 0) %>%
      st_as_sf(., coords = c('X', 'Y'), crs = 25833)
    randomp[[i]] <- p
  }

  # Get a dataframe out of the random points
  random_p <- do.call(rbind, randomp)
  
  random_p <- random_p %>% 
    mutate(FK_RegNr = id, season = NA, type = 'random_telem') %>% 
    dplyr::select(FK_RegNr, type, value)
  
  random_p$FK_RegNr <- as.numeric(random_p$FK_RegNr) 
  
  ###############################
  # Wrangle the CS observations #
  ###############################
  
  # Citizen science
  CS <- openxlsx::read.xlsx('real_analysis/dataset/cs/cs_mammals.xlsx', sheet = 1)
  CS$date <- paste(CS$year, CS$month, CS$day, sep ="-") %>% ymd() %>% as.Date() 
  
  moose_cs <- CS %>% 
    cutData(., type = 'season') %>% 
    filter(season == "summer (JJA)") %>% 
    filter(species == "Alces alces" & occurrenceStatus == "present") %>% 
    dplyr::select(species, month, decimalLongitude, decimalLatitude) %>% 
    st_as_sf(., coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>% 
    mutate(id = as.factor(1:nrow(.)), 
           type = 'cs',
           value = 1) %>% 
    st_transform(crs = 25833)
  # Unfortunately I do not have enough CS observations to afford filtering the
  # uncertainties
  
  # Take the Citizen Science observations inside the MCP buffer
  moose_cs <- st_intersection(moose_cs, mcp_buffer)
  
  # Select the variables I need
  moose_cs <- moose_cs %>% dplyr::select(id, value, type)
  
  ##########################################
  # Sample random citizen science absences #
  ##########################################
  moose_cs_sp <- as_Spatial(moose_cs)
  mcp_cs <- mcp(moose_cs_sp, percent = 100) %>% 
    st_as_sf()
  
  random_cs <- st_sample(mcp_buffer, size = ratio_random_cs*nrow(moose_cs)) %>% st_as_sf()
  
  coords_r <- st_coordinates(random_cs)
  
  random_cs <- random_cs %>% st_drop_geometry() %>% 
    mutate(X = coords_r[,1],
           Y = coords_r[,2]) %>%
    st_as_sf(., coords = c('X', 'Y'), crs = 25833) %>% 
    mutate(FK_RegNr = NA, season = NA, type = 'random_cs', value = 0) %>% 
    dplyr::select(FK_RegNr, type, value)
  
  #########################
  # Make the full dataset #
  #########################
  moose_full <- bind_rows(data_telem, random_p, moose_cs, random_cs)
  
  moose_full <- moose_full %>% 
    dplyr::select(FK_RegNr, value, type, geometry)
  
  ######################
  # Return the dataset #
  ######################
  return(moose_cs)
}

