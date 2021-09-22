######################################################################
##################### Make the Roe deer dataset ######################
######################################################################

make_reindeer_dataset <- function(buffer_for_cs, buffer_for_telem, thinning, ratio_random_cs){

######################################
# Wrangle the telemetry observations #
######################################

# Select the columns I want to extract // Data for summer only and for day only
telem_reindeer <- read_sf('real_analysis/dataset/reindeer_gps/gps_data_reindeer.shp') %>% 
  mutate(FK_RegNr = animal_id, 
         t = as.POSIXct(acquisitio), 
         value = 1, 
         type = 'telem_gps',
         X = st_coordinates(.[,1]),
         Y = st_coordinates(.[,2])) %>% 
  mutate(H = as.numeric(format(t, "%H"))) %>% 
  filter(H > 8 & H < 22) %>% 
  dplyr::select(FK_RegNr, geometry, t,  value, type, H, X, Y) %>% 
  st_transform(crs = 25833) %>% 
  filter(FK_RegNr != '3174') # Have a strange MCP ...

# Check if there is enough telemetry for each animal
id_telem <-  telem_reindeer %>% 
  st_drop_geometry() %>% 
  group_by(FK_RegNr) %>% 
  summarise(n = n()) %>% 
  filter(n >= 100) %>% 
  dplyr::select(FK_RegNr)

# One reindeer ID with only 1 observation -> filter it out
telem_reindeer <- telem_reindeer %>% 
  filter(FK_RegNr %in% id_telem$FK_RegNr)

###################################################################
# Thin the telemetry data to avoid spatial / temporal correlation #
###################################################################
vector_id <- unique(telem_reindeer$FK_RegNr)
reindeer_thinned <- list()

for(i in 1:length(vector_id)){
  dat_track <- telem_reindeer %>% dplyr::filter(FK_RegNr == paste(vector_id[i])) %>% 
    amt::make_track(.x = X, .y = Y, .t = t, key = FK_RegNr)
  dat_thinned <- amt::track_resample(dat_track, rate = thinning, tolerance = hours(1))
  reindeer_thinned[[i]] <- dat_thinned
}

telem_reindeer <- do.call(rbind, reindeer_thinned) %>%  
  dplyr::select(FK_RegNr = key, geometry) %>% 
  mutate(value = 1, type = 'telem_gps')

#############################################################
# Wrangle to get the minimum convex polygon for each animal #
#############################################################

reindeer_sp <- telem_reindeer %>% 
  dplyr::select(FK_RegNr, geometry) %>% 
  as_Spatial(.)

mcp <- mcp(reindeer_sp, percent = 99.5)
mcp_sf <- st_as_sf(mcp) %>% st_transform(., crs = 25833)
mcp_union <- st_union(mcp_sf) %>% st_buffer(buffer_for_telem)

# Crop the observations out of the polygon
data_telem <- st_intersection(telem_reindeer, mcp_union)

##################################
# Sample random telemetry points #
##################################

# Compute the number of points that I need for each animal
n_points_id <- telem_reindeer %>% 
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
  mcpb <- st_buffer(mcp, buffer_for_telem) 
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

cs_reindeer<- CS %>% 
  cutData(., type = 'season') %>% 
  filter(season == "summer (JJA)") %>% 
  filter(species == "Rangifer tarandus" & occurrenceStatus == "present") %>% 
  filter(recorder1 != 'Tord Bretten') %>% # Filter out Tord who is professional 
  dplyr::select(species, month, decimalLongitude, decimalLatitude, recorder1) %>% 
  st_as_sf(., coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>% 
  mutate(id = as.factor(1:nrow(.)), 
         type = 'cs',
         value = 1) %>% 
  st_transform(crs = 25833)
# Unfortunately I do not have enough CS observations to afford filtering the
# uncertainties

# Take the Citizen Science observations inside the reindeer areas
reindeer_areas <- read_sf('real_analysis/dataset/reindeer_areas/reindeer_areas.shp') %>% 
  st_transform(crs = 25833)

cs_reindeer <- st_intersection(cs_reindeer, mcp_union)

#################################################
# Random points for the Citizen Science dataset #
#################################################

random_cs <- st_sample(mcp_union, ratio_random_cs*nrow(cs_reindeer)) %>% st_as_sf() 

coords_r <- st_coordinates(random_cs)

random_cs <- random_cs %>% st_drop_geometry() %>% 
  mutate(X = coords_r[,1],
         Y = coords_r[,2]) %>%
  st_as_sf(., coords = c('X', 'Y'), crs = 25833) %>% 
  mutate(type = 'random_cs', value = 0) %>% 
  dplyr::select(type, value)

#########################
# Make the full dataset #
#########################
reindeer_full <- bind_rows(telem_reindeer, random_p, cs_reindeer, random_cs)

reindeer_full <- reindeer_full %>% 
  dplyr::select(FK_RegNr, value, type, geometry)

######################
# Return the dataset #
######################
return(reindeer_full)
}
