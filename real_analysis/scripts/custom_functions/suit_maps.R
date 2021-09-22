#########################################
# Function to make the suitability maps #
#########################################

make_suitability_map <- function(species, data, models, cell_size){

# Make the area for which to compute the grid
data_sp <- data %>% filter(type == 'telem_gps') %>% mutate(id = as.factor(FK_RegNr)) %>% 
  dplyr::select(id) %>% as_Spatial()

mcp_data <- mcp(data_sp) %>% st_as_sf() %>% st_union() %>% st_buffer(dist = 5000) %>% st_as_sf()

# Extract the coefficients
coef_model_telem <- models[[3]]$summary.fixed$mean
coef_model_naive <- models[[4]]$summary.fixed$mean
coef_model_corrected <- models[[5]]$summary.fixed$mean
coef_model_sp1 <- models[[6]]$summary.fixed$mean
coef_model_sp2 <- models[[7]]$summary.fixed$mean

# Make the grid
grid_data <- st_make_grid(mcp_data, cellsize = cell_size) %>% st_as_sf() %>% mutate(id = 1:nrow(.))

if(species == 'wreindeer'){
  random_points_data <- st_sample(grid_data, size = 100000) %>% st_as_sf() %>% get_covariate('wreindeer', ., list_raster_wreindeer, 1)}
else{
  random_points_data <- st_sample(grid_data, size = 100000) %>% st_as_sf() %>% get_covariate('other', ., list_raster_mr, 1)}

# Make the model matrix
intersect_data <- st_intersection(random_points_data, grid_data) %>% st_drop_geometry() %>% 
  group_by(id) %>% 
  summarise(slope = mean(slope),
            d_roads_log = mean(d_roads_log),
            d_urb_log = mean(d_urb_log),
            n_agr_log = mean(n_agr_log),
            n_forest_log = mean(n_forest_log),
            path_use_log = mean(path_use_log),
            alt_log = mean(alt_log)) %>% 
  full_join(., grid_data, by = 'id') %>% 
  drop_na()

Xmat <- model.matrix(~ d_roads_log + d_urb_log + path_use_log + slope + n_forest_log + n_agr_log + alt_log, data = intersect_data)

# Make the predictions
intersect_data <-intersect_data %>%  
  mutate(pred_telem = Xmat %*% coef_model_telem,
         pred_naive = Xmat %*% coef_model_naive,
         pred_corrected = Xmat %*% coef_model_corrected,
         pred_sp1 = Xmat %*% coef_model_sp1,
         pred_sp2 = Xmat %*% coef_model_sp2) %>% 
  st_as_sf()

# Plot the maps
a <- ggplot(intersect_data) + geom_sf(aes(fill = pred_telem), color = NA, show.legend = FALSE) + scale_fill_viridis(option = 'magma') + theme_map()
b <- ggplot(intersect_data) + geom_sf(aes(fill = pred_naive), color = NA, show.legend = FALSE) + scale_fill_viridis(option = 'magma') + theme_map()
c <- ggplot(intersect_data) + geom_sf(aes(fill = pred_corrected), color = NA, show.legend = FALSE) + scale_fill_viridis(option = 'magma') + theme_map()
d <- ggplot(intersect_data) + geom_sf(aes(fill = pred_sp1), color = NA, show.legend = FALSE) + scale_fill_viridis(option = 'magma') + theme_map()
e <- ggplot(intersect_data) + geom_sf(aes(fill = pred_sp2), color = NA, show.legend = FALSE) + scale_fill_viridis(option = 'magma') + theme_map()

#####################
# Objects to return #
#####################
maps_plot1 <- cowplot::plot_grid(a,b,c, nrow = 1)
cor_plot1 <- intersect_data %>% st_drop_geometry() %>% dplyr::select(pred_telem, pred_naive, pred_corrected) %>% drop_na() %>% cor(.) %>% ggcorrplot(., type = "lower", lab = TRUE)

maps_plot2 <- cowplot::plot_grid(a,b,c,d,e, nrow = 2)
cor_plot2 <- intersect_data %>% st_drop_geometry() %>% dplyr::select(pred_telem, pred_naive, pred_corrected, pred_sp1, pred_sp2) %>% drop_na() %>% cor(.) %>% ggcorrplot(., type = "lower", lab = TRUE)

###
return(list(maps_plot1,cor_plot1, maps_plot2, cor_plot2))
}



