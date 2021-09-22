######################################################################
##################### Make the Roe deer dataset ######################
######################################################################

make_roedeer_dataset <- function(buffer_for_cs, buffer_for_telem, thinning, ratio_random_cs){
  
  ######################################
  # Wrangle the telemetry observations #
  ######################################
  
  # Open the dataset and make the necessary transformations:
  # Change CRS
  # Make a "data" column readable by the amt package
  # Add a season index -> we will just need summer in our analysis
  telem_oslo <- read_sf('real_analysis/dataset/roedeer/Telem_RoeDeer_Oslo.shp') %>%
    dplyr::select(FK_RegNr = FK_RgNr, LOC_TIM, geometry) %>% 
    drop_na() %>% 
    st_as_sf() %>% 
    st_transform(., crs = 25833) %>% 
    dplyr::mutate(X = st_coordinates(.)[,1],
                  Y = st_coordinates(.)[,2],
                  t = as.POSIXct(LOC_TIM),
                  value = 1,
                  type = 'telem_gps') %>%   
    mutate(date = t) %>% 
    cutData(., type = 'season') %>% 
    dplyr::filter(season == "summer (JJA)") %>% 
    dplyr::select(FK_RegNr, geometry, X, Y, t, value, type)
  
  # We filter out the observations FK_Nr = X- as they are only tests
  # made by humans -> represent human locations
  
  # We also get rid of FK_Nr = 'R3063' as it messes with the computation
  # of the MCP
  telem_oslo <- telem_oslo  %>%
    dplyr::filter(!str_detect(FK_RegNr, "X")) %>% 
    dplyr::filter(FK_RegNr != 'R3063')
  
  # We get rid of FK_RegNr with less than 30 observations - to be sure they are
  # not other test GPS & to avoid risk of having dead animal still returning locations
  id_telem <-  telem_oslo %>% 
    st_drop_geometry() %>% 
    group_by(FK_RegNr) %>% 
    summarise(n = n()) %>% 
    filter(n >= 30) %>% 
    dplyr::select(FK_RegNr)
  
  # Filter out FK_Rg with less than 30 locations
  telem_oslo <- telem_oslo %>% 
    filter(FK_RegNr %in% id_telem$FK_RegNr)
  length(unique(telem_oslo$FK_RegNr))
  ###########################
  # Wrangle the VHF dataset #
  ###########################
  
  # We have 2 VHF dataset that we bind with the GPS data to increase the
  # coverage of the roe deer presence, thus increasing the number of
  # citizen science observations.
  
  # For both VHF dataset we:
  
  # - Give the dataset the good coordinates (25832) and transform it to crs = 3035
  # - Select observations on the ground & with an accuracy of less than 100m uncertainty
  # - We get rid of unreliable values (hours or minutes > 60 & days > 31)
  # - Correct the years adding "19"
  # - Create a "date" / t column that is understood by the amt package 
  # - select only the columns we are interested in
  # - Get only the summer observations
  
  ### Ostedata - 1st VHF dataset
  
  oste_data_select <- readxl::read_excel(path = "real_analysis/dataset/roedeer/VHF_data.xlsx", sheet = 1) %>% 
    drop_na(WGS84_X2) %>% 
    st_as_sf(., coords = c("WGS84_X2", "WGS84_Y2"), crs = 25832) %>% 
    st_transform(., crs = 25833) %>% 
    dplyr::filter(TYPE == "G") %>%
    dplyr::filter(ACC <= 100) %>% 
    dplyr::filter(H < 60) %>% 
    dplyr::filter(DD <= 31) %>% 
    dplyr::filter(DD >= 1) %>% 
    dplyr::filter(MIN < 60) %>% 
    mutate(year = paste0("19", YR)) %>% 
    mutate(year_num = as.numeric(year)) %>% 
    mutate(LOC_TIM = lubridate::make_datetime(year_num, MM, DD, H, MIN)) %>% 
    mutate(t = as.POSIXct(LOC_TIM)) %>% 
    mutate(FK_RegNr = NO, X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>% 
    dplyr::select(FK_RegNr, geometry, t, X, Y) %>% 
    dplyr::mutate(value = 1, type = 'telem_gps')
  
  # Get the summer observations
  oste_data_select <- oste_data_select %>% mutate(date = as.Date(t)) %>% 
    cutData(., type = 'season') %>% filter(season == 'summer (JJA)')
  
  ### Akerhus data - 2nd VHF dataset
  
  akerhus_data_select <- readxl::read_excel(path = "real_analysis/dataset/roedeer/VHF_data.xlsx", sheet = 2)  %>% 
    drop_na(WGS84X) %>% 
    st_as_sf(., coords = c("WGS84X", "WGS84Y"), crs = 25832) %>% 
    st_transform(., crs = 25833) %>% 
    dplyr::filter(HH < 60) %>%
    dplyr::filter(MN < 60) %>% 
    dplyr::filter(TYP == "G") %>% 
    dplyr::filter(POS <= 100) %>% 
    mutate(LOC_TIM = lubridate::make_datetime(YY, MO, DD, HH, MN)) %>% 
    mutate(t = as.POSIXct(LOC_TIM)) %>% 
    mutate(FK_RegNr = No., X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>% 
    dplyr::select(FK_RegNr, geometry, t, X, Y) %>% 
    dplyr::mutate(value = 1, type = 'telem_gps')
  
  # Select only the summer observations
  akerhus_data_select <- akerhus_data_select %>% mutate(date = as.Date(t)) %>% 
    cutData(., type = 'season') %>% filter(season == 'summer (JJA)')
  
  # Combine the two datasets
  vhf_bind <- bind_rows(oste_data_select, akerhus_data_select)
  
  # Filter out FK_RegNr with less than 20 locations
  id_vhf <- vhf_bind %>% 
    group_by(FK_RegNr) %>% 
    summarise(n = n()) %>% 
    filter(n >= 20) %>% 
    dplyr::select(FK_RegNr)
  
  # Select the relevant animals ID
  vhf_full <- vhf_bind %>% 
    filter(FK_RegNr %in% id_vhf$FK_RegNr)

  ###################################
  # Combine the 2 telemetry dataset #
  ###################################
  
  # We finally create a telemetry dataset combining both
  # the GPS telemetry data and the VHF data
  telem_oslo <- telem_oslo %>% 
    bind_rows(., vhf_full) %>% 
    st_as_sf() %>% 
    mutate(H = as.numeric(format(t, "%H"))) %>% 
    filter(H > 8 & H < 22)
  
  # We also get rid of observations which are outside the area of interest
  bbox <- read_sf("real_analysis/dataset/study_area/study_area.shp") %>% 
    st_transform(crs = 25833)
  telem_oslo <- st_crop(telem_oslo, bbox)
  
  ###################################################################
  # Thin the telemetry data to avoid spatial / temporal correlation #
  ###################################################################
  
  # To avoid any risk of spatial / temporal correlation we thin the dataset by taking 
  # a telemetry point every 5 hours. The time resolution of the current data is 2 hours
  
  vector_id <- unique(telem_oslo$FK_RegNr)
  telem_oslo_thinned <- list()
  
  # Run the subsampling algorithm
  for(i in 1:length(vector_id)){
    dat_track <- telem_oslo %>% dplyr::filter(FK_RegNr == paste(vector_id[i])) %>% 
      amt::make_track(.x = X, .y = Y, .t = t, all_cols = TRUE)
    dat_thinned <- amt::track_resample(dat_track, rate = thinning, tolerance = lubridate::hours(1))
    telem_oslo_thinned[[i]] <- dat_thinned
  }
  
  # Make it as a dataframe
  telem_oslo_thinned <- do.call(rbind, telem_oslo_thinned) %>%  
    dplyr::select(FK_RegNr, geometry, value) %>% 
    st_as_sf()
  
  #############################################################
  # Wrangle to get the minimum convex polygon for each animal #
  #############################################################
  
  # We transform the sf dataset into a sp dataset so it can be understood
  # by the mcp function
  for_mcp <- telem_oslo_thinned %>% 
    mutate(id = as.factor(FK_RegNr)) %>% 
    dplyr::select(id) %>% 
    as_Spatial()
  
  mcp <- mcp(for_mcp, percent = 97)

  # Merge the polygons and create a buffer around them so we cample
  # a bit more citizen science observations
  mcp_sf <- st_as_sf(mcp)
  mcp_union <- st_union(mcp_sf)
  mcp_buffer <- st_buffer(mcp_union, dist = buffer_for_cs) %>% 
    st_as_sf(crs = 25833)
  
  # Crop the observations out of the polygon
  data_telem <- st_intersection(telem_oslo_thinned, mcp_buffer)
  
  # Finally take the columns I need
  data_telem <- data_telem %>% 
    dplyr::select(FK_RegNr, geometry, value) %>% 
    mutate(type = "telem_gps")
  
  ######################################
  # Get random points telemetry points #
  ######################################
  
  # For each presence point we sample 2 absence points
  # Enough for the Infinitively Weighted Logistic Regression
  # (Muff et al. 2020)
  
  # Compute the number of points per animal ID
  n_points_id <- data_telem %>% 
    st_drop_geometry() %>%
    mutate(id = as.factor(FK_RegNr)) %>% 
    group_by(id) %>% 
    summarise(n = n()) %>% 
    mutate(to_sample = 2*n)
  
  mcp_sf <- mcp_sf %>% full_join(., n_points_id, by = 'id') %>% drop_na()
  
  # Make a list of mcps
  randomp <- list()
  
  # Run the sampling algorithm
  for(i in 1:nrow(mcp_sf)){
    mcp <- mcp_sf[i,]
    mcp <- st_buffer(mcp, buffer_for_telem)
    p <- mcp %>% st_sample(., mcp$to_sample) %>% st_as_sf()
    coords_p <- st_coordinates(p)
    p <- p %>% st_drop_geometry() %>% 
      mutate(X = coords_p[,1],
             Y = coords_p[,2],
             id = mcp$id,
             value = 0) %>%
      st_as_sf(., coords = c('X', 'Y'), crs = 25833)
    randomp[[i]] <- p
  }
  
  # Get a dataframe out of the random points, will be merged with the other data
  random_p <- do.call(rbind, randomp)
  
  random_p <- random_p %>% 
    mutate(FK_RegNr = id, season = NA, type = 'random_telem') %>% 
    dplyr::select(FK_RegNr, type, value)
  
  ###############################
  # Wrangle the CS observations #
  ###############################
  
  # Citizen science
  CS <- openxlsx::read.xlsx('real_analysis/dataset/cs/cs_mammals.xlsx', sheet = 1)
  CS$date <- paste(CS$year, CS$month, CS$day, sep ="-") %>% ymd() %>% as.Date() 
  
  cs_oslo <- CS %>% 
    cutData(., type = 'season') %>% 
    filter(season == "summer (JJA)") %>% 
    filter(species == "Capreolus capreolus" & occurrenceStatus == "present") %>% 
    dplyr::select(species, month, decimalLongitude, decimalLatitude) %>% 
    st_as_sf(., coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>% 
    mutate(id = as.factor(1:nrow(.)), 
           type = 'cs',
           value = 1) %>% 
    st_transform(crs = 25833)
  # Unfortunately I do not have enough CS observations to afford filtering the
  # uncertainties
  
  # Take the Citizen Science observations inside the MCP buffer
  cs_oslo <- cs_oslo %>% st_intersection(., mcp_buffer)
  
  #################################################
  # Random points for the Citizen Science dataset #
  #################################################
  rd_cs_sp <- as_Spatial(cs_oslo)
  mcp_cs <- mcp(rd_cs_sp, percent = 100) %>% 
    st_as_sf()
  
  random_cs <- st_sample(mcp_buffer, size = ratio_random_cs*nrow(cs_oslo)) %>% st_as_sf()
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
  
  # Finally make the full dataset ...
  data <- bind_rows(cs_oslo, data_telem, random_p, random_cs) %>% 
    dplyr::select(FK_RegNr, type, value, geometry)
  
  ######################
  # Return the dataset #
  ######################
  return(data)
}

