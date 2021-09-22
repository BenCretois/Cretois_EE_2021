# R code for the paper Cretois et al. 'Identifying and correcting spatial bias in opportunistic citizen science data for wild ungulates in Norway'

---

This repository hosts functions . 

**Note that to run the "real analysis"** the raster used to extract the covariates, the telemetry and Citizen Science dataset were too heavy to be shared on GitHub. You can find all this data at https://zenodo.org/record/4590153#.YGGdMq8zZaQ. **Place the downloaded folders in data to smoothly run the scripts**

This work was supported by a PhD grant funded by the Norwegian University of Science and Technology and the Research Council of Norway.

---

## Simulation study

[simul_analysis](./simulation_study/simul_analysis.R): Script running the entire simulation.

In the script [simul_functions](./simulation_study/simul_functions.R) are the functions used to run the simulation:

* simul_env: function simulating all the environmental variables
* simul_data: function simulating both the species' presence and the citizen science observations
* fit_models: function fitting a RSF using the species presence wit randomly sampled availability points, a RSF using citizen science observations with randomly sampled availability points and a RSF using citizen science observations with corrected availability points
* boot_rmse: a function computing bootstraping the difference in RMSE between the naive CS model and the corrected CS model

More details concerning the simulation can be found in the Annexes of the paper.

[bootstrap_data_gen.R](.simulation_study/bootstrap_data_gen.R) aims at running a full analysis to allow bootstrap of the data generation
process as well as running the model.


## Analysis with real data


To run the whole analysis without going through the individual functions you can use the [run_all](./real_study/run_all.R) functions.

### Data wrangling

Functions below aim to create a simplified and standardized dataset between species.

[make_moose_dataset](./real_study/custom_functions/make_moose_dataset.R): Open the moose GPS telemetry data, the citizen science dataset and create a simplified dataset that contains only summer observations. Buffer around moose minimum convex polygon can be specified to sample more or less citizen science observations. Thinning rate can also be specified to include more or less GPS telemetry locations.

[make_roedeer_dataset](./real_study/custom_functions/make_roedeer_dataset.R): Open the roe deer GPS telemetry data and VHF data, the citizen science dataset and create a simplified dataset that contains only summer observations. Buffer around moose minimum convex polygon can be specified to sample more or less citizen science observations. Thinning rate can also be specified to include more or less GPS telemetry locations.

[make_reindeer_dataset](./real_study/custom_functions/make_moose_dataset.R): Open the wild reindeer GPS telemetry data, the citizen science dataset and create a simplified dataset that contains only summer observations. Buffer around moose minimum convex polygon can be specified to sample more or less citizen science observations. Thinning rate can also be specified to include more or less GPS telemetry locations.

[get_covariate](./real_study/custom_functions/get_covariate.R): Extract the raster value for each points in the dataset (GPS telemetry, citizen science and availability).

[covariate_plot](./real_study/custom_functions/covariate_plot.R): Function to plot the distribution of the covariate value for GPS telemetry, citizen science data and availability

### Sample a corrected availability 

[get_corrected_availability](./real_study/custom_functions/get_corrected_availability.R): Internally, the function fits an observer model using telemetry data (coded as 0) and citizen science observations (coded as 1). Once the observer model is fitted, the function draw a map of probability of having a citizen science observation at a specific grid cell. The function then sample availability points weighted by this probability (In [run_all](./real_study/run_all.R) we use *ratio_availability* = 2, meaning that we sample 2*the number of citizen science observation) 

[get_availability_cross_species](./real_study/custom_functions/get_availability_cross_species.R): The function takes the observer model of another species as an argument and use it to sample the corrected availability points (procedure similar to the function [get_corrected_availability](./real_study/custom_functions/get_corrected_availability.R))

[fit_models](./real_study/custom_functions/fit_models.R): The function fits three RSF, 
* A model using the species' GPS telemetry location with random availability across the species' range
* A model using citizen science observations and random availability across the species' range 
* A model using citizen science observation with corrected availability using the target species' observer model
* Two model using citizen science observation with corrected availability using the other two species' observer model

### Displaying the results

[make_suitability_map](./real_study/custom_functions/suit_maps.R): The function uses the coefficient of the fitted RSF to make habitat preference maps.



