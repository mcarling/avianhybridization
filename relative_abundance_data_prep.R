### In this script, we prepare datasets for modeling the relative abundances of 
### two parental species and hybrids. First, we calculate habitat, elevation,
### and bioclim variables for all checklists in the eBird dataset. We then 
### generate a prediction surface of these variables for North America.

# loading required packages
library(sf)
library(raster)
library(MODIS)
library(exactextractr)
library(tidyverse)
library(lubridate)

# additional packages needed to add bioclim data
library(sp)
library(maptools)
library(rgdal)
library(dismo)

# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

##### loading filtered and zero-filled eBird data -----------------------------

setwd("/Volumes/commons/CarlingLab/eBird Data/Data for looking at relative abundance")

## for this project, we filtered down to all checklists submitted in Canada, US, and Mexico from June 1st to July 16th
## Note: because the eBird data for our three taxa of interest have already been zero-filled, the locality data
## for each should be identical. Therefore, we can generate habitat, elevation, and bioclim data for checklist 
## locations for only one dataset. Here, I'll do it for the Indigo Bunting data
ebird <- read.csv("indigo_bunting_pred_mx_half_july.csv", header = TRUE)

###### preparing MODIS landcover data ------------------------------------------
### creating a landcover classification around each point of the eBird dataset
# load the landcover data

## making sure to filter down ne_land shapefile to only land that 
## we want to plot in our final maps
map_proj <- st_crs(102003) # the projection we'll use for plotting
ne_land <- read_sf("data_eBird_best_practices/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj)

# crop down to an area that comprises the (more efficient to generate
# prediciton surface data for this subset than the entire continent)
ne_land_trim_final <- st_crop(ne_land, c(ymin=-930000, ymax=1500000, xmin=-2500000, xmax=2150000))

# now converting to native modis projection
ne_land_modis <- ne_land_trim_final %>% 
  # project to the native modis projection
  st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))
plot(ne_land_modis)


tiles <- getTile(ne_land_modis)
tiles@tile

# earliest year of ebird data
begin_year <- "2010.01.01"
# end date for ebird data
end_year <-  "2018.12.31"
# download tiles and combine into a single raster for each year
tifs <- runGdal(product = "MCD12Q1", collection = "006", SDSstring = "01",
                extent = ne_land_modis %>% st_buffer(dist = 10000), # probably don't need to add the buffer if downloading data for all of North America, right?
                begin = begin_year, end = end_year, 
                outDirPath = "modis_landcover2", job = "modis",
                MODISserverOrder = "LPDAAC") %>% 
  pluck("MCD12Q1.006") %>% 
  unlist()

# rename tifs to have more descriptive names
new_names <- format(as.Date(names(tifs)), "%Y") %>% 
  sprintf("modis_mcd12q1_umd_%s.tif", .) %>% 
  file.path(dirname(tifs), .)
file.rename(tifs, new_names)


# load the landcover data
landcover <- list.files("modis_landcover2/modis", "^modis_mcd12q1_umd", 
                        full.names = TRUE) %>% 
  stack()
# label layers with year
landcover <- names(landcover) %>% 
  str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover, .)
landcover


##### combining MODIS and eBird data ------------------------------------------

#converting observation date to date, not a factor
ebird$observation_date <- as_date(ebird$observation_date)

max_lc_year <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  as.integer() %>% 
  max()

# creating a 2.5 buffer around each point
neighborhood_radius <- 5 * ceiling(max(res(landcover))) / 2
ebird_buff <- ebird %>%
  distinct(year = format(observation_date, "%Y"),
           locality_id, latitude, longitude) %>% 
  # for 2019 and 2020, use 2018 landcover data
  mutate(year_lc = if_else(as.integer(year) > max_lc_year,
                           as.character(max_lc_year), year),
         year_lc = paste0("y", year_lc)) %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to modis projection
  st_transform(crs = projection(landcover)) %>% 
  # buffer to create neighborhood around each point
  st_buffer(dist = neighborhood_radius) %>% 
  # nest by year
  nest(data = c(year, locality_id, geometry))

# function to summarize landcover data for all checklists in a given year
calculate_pland <- function(yr, regions, lc) {
  locs <- st_set_geometry(regions, NULL)
  exact_extract(lc[[yr]], regions, progress = FALSE) %>% 
    map(~ count(., landcover = value)) %>% 
    tibble(locs, data = .) %>% 
    unnest(data)
}

# iterate over all years extracting landcover for all checklists in each
ebird_lc_extract <- ebird_buff %>%
  mutate(pland = map2(year_lc, data, calculate_pland, lc = landcover)) %>% 
  select(pland) %>% 
  unnest(cols = pland)

## calculating PLAND for each observation location
ebird_pland <- ebird_lc_extract %>%
  # calculate proporiton
  group_by(locality_id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

# convert names to be more descriptive
lc_names <- tibble(landcover = 0:15,
                   lc_name = c("pland_00_water", 
                               "pland_01_evergreen_needleleaf", 
                               "pland_02_evergreen_broadleaf", 
                               "pland_03_deciduous_needleleaf", 
                               "pland_04_deciduous_broadleaf", 
                               "pland_05_mixed_forest",
                               "pland_06_closed_shrubland", 
                               "pland_07_open_shrubland", 
                               "pland_08_woody_savanna", 
                               "pland_09_savanna", 
                               "pland_10_grassland", 
                               "pland_11_wetland", 
                               "pland_12_cropland", 
                               "pland_13_urban", 
                               "pland_14_mosiac", 
                               "pland_15_barren"))

ebird_pland <- ebird_pland %>%
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)

ebird_pland <- ebird_pland %>%
  pivot_wider(names_from = lc_name, 
              values_from = pland, 
              values_fill = list(pland = 0))

write_csv(ebird_pland, "ebird_pland_location_final,year.csv", na = "") # for Lazuli
ebird_pland<-read.csv("ebird_pland_location_final,year.csv")

##### making MODIS prediction surface for models ------------------------------

agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
  aggregate(agg_factor) 
r <- ne_land_modis %>%  #doing it for all of North America, not just BCR27
  st_transform(crs = projection(r)) %>%
  rasterize(r, field = 1) %>%
  # remove any empty cells at edges
  trim()

r <- writeRaster(r, filename = "prediction-surface_final.tif", overwrite = TRUE)
r <- raster("prediction-surface_final.tif") 

# get cell centers and create neighborhoods
r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
  st_as_sf() %>% 
  transmute(id = row_number())
r_cells <- st_buffer(r_centers, dist = neighborhood_radius)

# extract landcover values within neighborhoods, only needed most recent year
# note: this step can take a while
lc_extract_pred <- landcover[[paste0("y", max_lc_year)]] %>% 
  exact_extract(r_cells, progress = FALSE) %>% 
  map(~ count(., landcover = value)) %>% 
  tibble(id = r_cells$id, data = .) %>% 
  unnest(data)

# calculate the percent for each landcover class
pland_pred <- lc_extract_pred %>% 
  count(id, landcover) %>% 
  group_by(id) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

# convert names to be more descriptive
pland_pred <- pland_pred %>% 
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)

# tranform to wide format, filling in implicit missing values with 0s
pland_pred <- pland_pred %>% 
  pivot_wider(names_from = lc_name, 
              values_from = pland, 
              values_fill = list(pland = 0)) %>% 
  mutate(year = max_lc_year) %>% 
  select(id, year, everything())

# join in coordinates
pland_coords <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred, by = "id")

## adding elevation values
## see directions on best practices for using eBird guide for how 
## to download EarthEnv elevation data
elev <- raster("elevation_1KMmd_GMTEDmd.tif")
# crop, per Strimas-Mackey et al, buffering cropped ne_land shapefile by 10 km
elev <- ne_land_modis %>% 
  st_buffer(dist = 10000) %>% 
  st_transform(crs = projection(elev)) %>% 
  crop(elev, .) %>% 
  projectRaster(crs = projection(landcover))

## Note: for adding elevation values for eBird dataset,
## it might be best to download elevation data for all North America (optional)
# ne_land <- read_sf("data_eBird_best_practices/gis-data.gpkg", "ne_land") %>% 
#   # project to the native modis projection
#   st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
#                            "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))
# elev <- ne_land %>% 
#   st_buffer(dist = 10000) %>% 
#   st_transform(crs = projection(elev)) %>% 
#   crop(elev, .) %>% 
#   projectRaster(crs = projection(landcover))

## buffer each checklist location
ebird_buff_noyear <- ebird %>%
  distinct(locality_id, latitude, longitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(elev)) %>% 
  st_buffer(dist = neighborhood_radius)

# extract elevation values and calculate median and sd
ebird_locs <- st_set_geometry(ebird_buff_noyear, NULL) %>% # for Lazulis
  mutate(id = row_number())
ebird_elev_checklists <- exact_extract(elev, ebird_buff_noyear, progress = FALSE) %>% 
  map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                   elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
  # join to lookup table to get locality_id
  bind_cols(ebird_locs, .)

# calculating elevation values for prediction surface
# extract and calculate median and sd
## this step can take a really, really long time
elev_pred <- exact_extract(elev, r_cells, progress = FALSE) %>% 
  map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE))) %>% #, # not doing sd to hopefully expedite this
                 #  elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
  # join to lookup table to get locality_id
  bind_cols(st_drop_geometry(r_cells), .)

### combining elevation and landcover covariates in the same datasets
# checklist covariates
ebird_pland_elev_checklist <- inner_join(ebird_pland, ebird_elev_checklists, by = "locality_id")
write_csv(ebird_pland_elev_checklist, "ebird_pland-elev_location-year_final.csv")

# prediction surface covariates
pland_elev_pred <- inner_join(pland_coords, elev_pred, by = "id")
write_csv(pland_elev_pred, "pland-elev_prediction-surface_trim_final.csv")

# adding bioclim data to eBird checklists and prediction surface ----------

# organizing workspace
dir.create(path = "data")
dir.create(path = "output")

# accessing worldclim data
bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5, # can change resolution to 2.5, 5, or 10 minutes of a degree
                        path = "data/")

# creating a dataframe with the full name of each bioclim variable
bioclim_names <- tibble(bioclim_var = paste("bio", 1:19, sep = ""),
                        bioclim_name = c("bio1_annual_mean_temp",
                                         "bio2_mean_diurnal_range",
                                         "bio3_isothermality",
                                         "bio4_temperature_seasonality",
                                         "bio5_max_temp_warmest_period",
                                         "bio6_min_temp_coldest_period",
                                         "bio7_temp_annual_range",
                                         "bio8_mean_temp_wettest_quarter",
                                         "bio9_mean_temp_driest_quarter",
                                         "bio10_mean_temp_warmest_quarter",
                                         "bio11_mean_temp_coldest_quarter",
                                         "bio12_annual_precip",
                                         "bio13_precip_wettest_period",
                                         "bio14_precip_driest_period",
                                         "bio15_precip_seasonality",
                                         "bio16_precip_wettest_quarter",
                                         "bio17_precip_driest_quarter",
                                         "bio18_precip_warmest_quarter",
                                         "bio19_precip_coldest_quarter"))


### adding bioclim data to eBird checklists

ebird_sf <- ebird %>% # for hybrids
  # filtering down to unique localities, don't care about year for this, I don't think
  # we also don't need to create a buffer around each point, all we need is the value for each point
  distinct(locality_id, latitude, longitude) %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% # have to change crs? no, you shouldn't need to
  # transform to bioclim projection
  st_transform(crs = projection(bioclim.data))

# ## buffer each checklist location
# ebird_buff_noyear <- ebird %>%
#   distinct(locality_id, latitude, longitude) %>% 
#   st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
#   st_transform(crs = projection(bioclim.data)) %>% 
#   st_buffer(dist = neighborhood_radius)


# need to pull out each variable from the bioclim dataset to add to lazuli sf dataset!!!!!!!
ebird_locs <- st_set_geometry(ebird_sf, NULL) %>% # for Lazulis
  mutate(id = row_number())

relevant_vars <- c(10,12,18)

bioclim_values = vector('list', length(relevant_vars))

for (i in relevant_vars){
# for now not creating a buffer
  bioclim_checklists <- raster::extract(bioclim.data[[i]], ebird_sf)
  bioclim_values[[i]] <- tibble(locality_id = ebird_locs$locality_id, bioclim_var = names(bioclim.data[[i]]), bioclim_value = bioclim_checklists)
  
}
bioclim_values_all <- do.call(rbind, bioclim_values)

# renaming bioclim variable names

# combining with new names
ebird_bioclim <- bioclim_values_all %>%
  inner_join(bioclim_names, by = "bioclim_var") %>% 
  arrange(bioclim_var) %>% 
  select(-bioclim_var)

ebird_bioclim <- ebird_bioclim %>% # for hybrid
  pivot_wider(names_from = bioclim_name, 
              values_from = bioclim_value, 
              values_fill = list(bioclim_value = NA))

ebird_pland_elev_bioclim_checklist <- inner_join(ebird_pland_elev_checklist, ebird_bioclim, by = "locality_id")
write_csv(ebird_pland_elev_bioclim_checklist, "ebird_pland-elev-bioclim_location-year_final.csv")


### adding bioclim data to prediction surface

relevant_vars <- c(10,12,18) # only looking at bioclim variables 
# that we know to relevant from Carling & Thomassen (2011)

# bioclim_values_pred = vector('list', nlayers(bioclim.data))
bioclim_values_pred = vector('list', length(relevant_vars))

for (i in relevant_vars){
  bioclim <- ne_land_modis %>%  
    st_buffer(dist = 10000) %>% # same as what we did with elevation
    st_transform(crs = projection(bioclim.data[[i]])) %>% # making sure projection is the same as the rest of the prediction surface
    crop(bioclim.data[[i]], .) %>% 
    projectRaster(crs = projection(landcover))
  
  bioclim_values_pred[[i]] <- exact_extract(bioclim, r_cells, progress = FALSE) %>% 
    map_dfr(~ tibble(bioclim_median = mean(.$value, na.rm = TRUE), # not doing sd to hopefully expedite this
                     bioclim_var = names(bioclim.data[[i]]))) %>% 
    # join to lookup table to get locality_id
    bind_cols(st_drop_geometry(r_cells), .)
}

# making sure the above for loop worked
head(bioclim_values_pred[[10]])
head(bioclim_values_pred[[12]])
head(bioclim_values_pred[[18]])

# combining all into a single dataframe
bioclim_values_pred_all <- do.call(rbind, bioclim_values_pred)

unique(bioclim_values_pred_all$bioclim_var) # making sure the above 
# steps worked, you should have the names of the bioclime variables 
# you're looking at

# combining with new names and spreading
bioclim_values_pred_all <- bioclim_values_pred_all %>%
  inner_join(bioclim_names, by = "bioclim_var") %>% 
  arrange(bioclim_var) %>% 
  select(-bioclim_var)
head(bioclim_values_pred_all)

bioclim_values_pred_all <- bioclim_values_pred_all %>%
  pivot_wider(names_from = bioclim_name, 
              values_from = bioclim_median, 
              values_fill = list(bioclim_value = NA))
head(bioclim_values_pred_all)

# prediction surface covariates
pland_elev_bioclim_pred <- inner_join(pland_elev_pred, bioclim_values_pred_all, by = "id")
write_csv(pland_elev_bioclim_pred, "pland-elev-bioclim_prediction-surface_final.csv")





##### Modeling relative abundance ---------------------------------------------

# set random number seed to insure fully repeatable results
set.seed(1)

# loading eBird data
ebird <- read_csv("lazuli_bunting_pred.csv") %>% 
  mutate(protocol_type = factor(protocol_type, 
                                levels = c("Stationary" , "Traveling"))) %>%
  # remove observations with no count
  filter(!is.na(observation_count)) %>%
  select(!X1)


# modis habitat covariates
habitat <- read_csv("lazuli_pland-elev_location-year.csv") %>% 
  mutate(year = as.integer(year))

# combine ebird and habitat data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))


# prediction surface
pred_surface <- read_csv("pland-elev_prediction-surface.csv")
# latest year of landcover data
max_lc_year <- pred_surface$year[1]
r <- raster("prediction-surface8.tif")


### skipping loading GIS data for making maps

## spatiotemporal subsampling
# generate hexagonal grid with ~ 5 km betweeen cells
dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id and week number for each checklist
checklist_cell <- ebird_habitat %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         week = week(observation_date))
# sample one checklist per grid cell per week
# sample detection/non-detection independently 
ebird_ss <- checklist_cell %>% 
  group_by(species_observed, year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell, -week)

## test train split
# including all habitat covariates here
hab_covs <- c("pland_00_water", 
              "pland_01_evergreen_needleleaf", 
              "pland_02_evergreen_broadleaf", 
              "pland_03_deciduous_needleleaf", 
              "pland_04_deciduous_broadleaf", 
              "pland_05_mixed_forest",
              "pland_06_closed_shrubland", 
              "pland_07_open_shrubland", 
              "pland_08_woody_savanna", 
              "pland_09_savanna", 
              "pland_10_grassland", 
              "pland_11_wetland", 
              "pland_12_cropland", 
              "pland_13_urban", 
              "pland_14_mosiac", 
              "pland_15_barren")
ebird_split <- ebird_ss %>% 
  # select only the columns to be used in the model
  select(observation_count,
         # effort covariates
         day_of_year, time_observations_started, duration_minutes,
         effort_distance_km, number_observers, protocol_type,
         # habitat covariates
         hab_covs,
         # also adding latitude, longitude, and elevation
         latitude, longitude, elevation_median)

# #scaling relevant predictor variables (is it right to do this all at the beginning, or do I need to scale training dataset first and then apply same sd to test dataset?)
# #ebird_split_scaled <- scale(ebird_split[, !ebird_split$protocol_type])
# ebird_split_scaled <- ebird_split %>%
#   mutate_at(c(2:6,8:26), funs(c(scale(.))))


# split 80/20
ebird_split <- ebird_split %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
map_int(ebird_split, nrow)


## exploratory data analysis
p <- par(mfrow = c(1, 2))
# counts with zeros
hist(ebird_ss$observation_count, main = "Histogram of counts", 
     xlab = "Observed count")
# counts without zeros
pos_counts <- keep(ebird_ss$observation_count, ~ . > 0)
hist(pos_counts, main = "Histogram of counts > 0", 
     xlab = "Observed non-zero count")
par(p)


## abundance models
# gam parameters
# degrees of freedom for smoothing
k <- 5
# degrees of freedom for cyclic time of day smooth
k_time <- 7 

# continuous predictors
# hold out time to treat seperately since it's cyclic
continuous_covs <- ebird_split$train %>% 
  select(-observation_count, -protocol_type, -time_observations_started) %>% 
  names()

# create model formula for predictors
gam_formula_rhs <- str_glue("s({var}, k = {k})", 
                            var = continuous_covs, k = k) %>% 
  str_flatten(collapse = " + ") %>% 
  str_glue(" ~ ", .,
           " + protocol_type + ",
           "s(time_observations_started, bs = \"cc\", k = {k})", 
           k = k_time) %>% 
  as.formula()

# model formula including response
gam_formula <- update.formula(observation_count ~ ., gam_formula_rhs)
gam_formula
#> observation_count ~ s(day_of_year, k = 5) + s(duration_minutes, 
#>     k = 5) + s(effort_distance_km, k = 5) + s(number_observers, 
#>     k = 5) + s(pland_04_deciduous_broadleaf, k = 5) + s(pland_05_mixed_forest, 
#>     k = 5) + s(pland_12_cropland, k = 5) + s(pland_13_urban, 
#>     k = 5) + protocol_type + s(time_observations_started, bs = "cc", 
#>     k = 7)
# try fewer knots for longitude?

### try scaling predictors!! actually maybe just scale latitude, longitude, and elevation
# gam_formula <- observation_count ~ s((scale(day_of_year)), k = 5) + 
#   s((scale(duration_minutes)), k = 5) + 
#   s(scale(effort_distance_km), k = 5) + 
#   s((scale(number_observers)), k = 5) + 
#   s((scale(pland_01_evergreen_needleleaf)), k = 5) + 
#   s((scale(pland_04_deciduous_broadleaf)), k = 5) + 
#   s((scale(pland_06_closed_shrubland)), k = 5) + 
#   s((scale(pland_07_open_shrubland)), k = 5) + 
#   s((scale(pland_08_woody_savanna)), k = 5) + 
#   s((scale(pland_09_savanna)), k = 5) +
#   s((scale(pland_12_cropland)), k = 5) + 
#   s((scale(pland_13_urban)), k = 5) +
#   s((scale(latitude)), k = 5) +
#   s((scale(longitude)), k = 5) +
#   s((scale(elevation_median)), k = 5) +
#   protocol_type + 
#   s((scale(time_observations_started)), bs = "cc", k = 7)

# specifying gam formula manually, scaling latitude, longitude, and elevation
gam_formula <- observation_count ~ s(day_of_year, k = 5) + s(duration_minutes, k = 5) + 
  s(effort_distance_km, k = 5) + s(number_observers, k = 5) + 
  s(pland_00_water, k = 5) + s(pland_01_evergreen_needleleaf, k = 5) + 
  s(pland_02_evergreen_broadleaf, k = 5) + s(pland_03_deciduous_needleleaf, k = 5) + 
  s(pland_04_deciduous_broadleaf, k = 5) + s(pland_05_mixed_forest, k = 5) + 
  s(pland_06_closed_shrubland, k = 5) + s(pland_07_open_shrubland, k = 5) + 
  s(pland_08_woody_savanna, k = 5) + s(pland_09_savanna, k = 5) + 
  s(pland_10_grassland, k = 5) + s(pland_11_wetland, k = 5) + 
  s(pland_12_cropland, k = 5) + s(pland_13_urban, k = 5) + 
  s(pland_14_mosiac, k = 5) + s(pland_15_barren, k = 5) + 
  s((scale(latitude)), k = 5) + s((scale(longitude)), k = 5) + s((scale(elevation_median)), k = 5) + 
  protocol_type + s(time_observations_started, bs = "cc", k = 7)

# trying to do the model with different error distributions
# explicitly specify where the knots should occur for time_observations_started
# this ensures that the cyclic spline joins the variable at midnight
# this won't happen by default if there are no data near midnight
time_knots <- list(time_observations_started = seq(0, 24, length.out = k_time))

# zero-inflated poisson
m_ziplss <- gam(list(gam_formula,      # count model
                     gam_formula[-2]), # presence model
                data = ebird_split$train, 
                family = "ziplss", 
                knots = time_knots)

# negative binomial
m_nb <- gam(gam_formula,
            data = ebird_split$train, 
            family = "nb",
            knots = time_knots)

# tweedie distribution
m_tw <- gam(gam_formula,
            data = ebird_split$train, 
            family = "tw",
            knots = time_knots)


### skipping model assesment for now


### examining covariate effects
# ggplot function
plot_gam <- function(m, title = NULL, ziplss = c("presence", "abundance")) {
  # capture plot
  tmp <- tempfile()
  png(tmp)
  p <- plot(m, pages = 1)
  dev.off()
  unlink(tmp)
  
  # drop addition models in ziplss
  if (m$family$family == "ziplss") {
    is_presence <- map_lgl(p, ~ str_detect(.$ylab, "^s\\.1"))
    if (ziplss == "presence") {
      p <- p[is_presence]  
    } else {
      p <- p[!is_presence]
    }
  }
  
  # extract data
  p_df <- map_df(p, ~ tibble(cov = rep(.$xlab, length(.$x)),
                             x = .$x, fit = .$fit, se = .$se))
  
  # plot
  g <- ggplot(p_df) +
    aes(x = x, y = fit,
        ymin = fit - se, ymax = fit + se) +
    geom_ribbon(fill = "grey80") +
    geom_line(col = "blue") +
    facet_wrap(~ cov, scales = "free_x") +
    labs(x = NULL,
         y = "Smooth function",
         title = title)
  print(g)
  invisible(p_df)
}

plot_gam(m_nb, title = "Negative Binomial GAM")


# defining prediction model
pred_model <- m_nb


### prediction
# create a dataframe of covariates with a range of start times
#seq_tod <- seq(0, 24, length.out = 300)
seq_tod <- seq(-2.737565, 2.887304, length.out = 300)
tod_df <-ebird_split_scaled$train %>% 
  # find average pland habitat covariates
  select(latitude, longitude, elevation_median, (starts_with("pland"))) %>% 
  summarize_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  # use standard checklist
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling") %>% 
  cbind(time_observations_started = seq_tod)

# predict at different start times
pred_tod <- predict(pred_model, newdata = tod_df, 
                    type = "link", 
                    se.fit = TRUE) %>% 
  as_tibble() %>% 
  # calculate backtransformed confidence limits
  transmute(time_observations_started = seq_tod,
            pred = pred_model$family$linkinv(fit),
            pred_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
            pred_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit))

# find optimal time of day
t_peak <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]

# plot the partial dependence plot
ggplot(pred_tod) +
  aes(x = time_observations_started, y = pred,
      ymin = pred_lcl, ymax = pred_ucl) +
  geom_ribbon(fill = "grey80", alpha = 0.5) +
  geom_line() +
  geom_vline(xintercept = t_peak, color = "blue", linetype = "dashed") +
  labs(x = "Hours since midnight",
       y = "Predicted relative abundance",
       title = "Effect of observation start time on Lazuli Bunting reporting",
       subtitle = "Peak detectability shown as dashed blue line")


# add effort covariates to prediction surface
pred_surface_eff <- pred_surface %>% 
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
         time_observations_started = t_peak,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling")

# predict
pred <- predict(pred_model, newdata = pred_surface_eff, 
                type = "link", 
                se.fit = TRUE) %>% 
  as_tibble() %>% 
  # calculate confidence limits and back transform
  transmute(abd = pred_model$family$linkinv(fit),
            abd_se = pred_model$family$linkinv(se.fit),
            abd_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
            abd_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit)) %>%
  # add to prediction surface
  bind_cols(pred_surface_eff, .) %>% 
  select(latitude, longitude, abd, abd_se, abd_lcl, abd_ucl)



### converting object to spatial dataframe
r_pred <- pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
r_pred <- r_pred[[-1]]

# save the rasters
# tif_dir <- "output"
# if (!dir.exists(tif_dir)) {
#   dir.create(tif_dir)
# }
# writeRaster(r_pred[["abd"]], 
#             filename = file.path(tif_dir, "abundance-model_abd_woothr.tif"),
#             overwrite = TRUE)
# writeRaster(r_pred[["abd_se"]], 
#             filename = file.path(tif_dir, "abundance-model_se_woothr.tif"), 
#             overwrite = TRUE)





#### making the relative abundance map!!!!

# any expected abundances below this threshold are set to zero
zero_threshold <- 0.05

# project predictions
map_proj <- st_crs(102003)
r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")



# load gis data for making maps
map_proj <- st_crs(102003)
ne_land <- read_sf("~/Downloads/data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

ne_land_trim2 <- ne_land_trim %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf("~/Downloads/data/gis-data.gpkg", "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("~/Downloads/data/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf("~/Downloads/data/gis-data.gpkg", "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()



par(mfrow = c(2, 1))
for (nm in names(r_pred)) {
  r_plot <- r_pred_proj[[nm]]
  
  par(mar = c(3.5, 0.25, 0.25, 0.25))
  # set up plot area
  plot(ne_land_trim2, col = NA, border = NA)
  plot(ne_land_trim2, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  
  # modified plasma palette
  plasma_rev <- rev(plasma(25, end = 0.9))
  gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
  pal <- c(gray_int(4)[2], plasma_rev)
  
  # abundance vs. se
  if (nm == "abd") {
    title <- "Lazuli Bunting Relative Abundance"
    # set very low values to zero
    r_plot[r_plot <= zero_threshold] <- NA
    # log transform
    r_plot <- log10(r_plot)
    # breaks and legend
    mx <- ceiling(100 * cellStats(r_plot, max)) / 100
    mn <- floor(100 * cellStats(r_plot, min)) / 100
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- sort(c(-2:2, mn, mx))
    lbls <- round(10^lbl_brks, 2)
  } else {
    title <- "Lazuli Bunting Abundance Uncertainty (SE)"
    # breaks and legend
    mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
    mn <- floor(1000 * cellStats(r_plot, min)) / 1000
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- seq(mn, mx, length.out = 5)
    lbls <- round(lbl_brks, 2)
  }
  
  # abundance
  plot(r_plot, 
       col = pal, breaks = brks, 
       maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE)
  
  # borders
  #plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()
  
  # legend
  par(new = TRUE, mar = c(0, 0, 0, 0))
  image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
             smallplot = c(0.25, 0.75, 0.06, 0.09),
             horizontal = TRUE,
             axis.args = list(at = lbl_brks, 
                              labels = lbls,
                              fg = "black", col.axis = "black",
                              cex.axis = 0.75, lwd.ticks = 0.5,
                              padj = -1.5),
             legend.args = list(text = title,
                                side = 3, col = "black",
                                cex = 1, line = 0))
}







par(mfrow = c(2, 1))
for (nm in names(r_pred)) {
  r_plot <- r_pred_proj[[nm]]
  
  par(mar = c(3.5, 0.25, 0.25, 0.25))
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  
  # modified plasma palette
  plasma_rev <- rev(plasma(25, end = 0.9))
  gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
  pal <- c(gray_int(4)[2], plasma_rev)
  
  # abundance vs. se
  if (nm == "abd") {
    title <- "Wood Thrush Relative Abundance"
    # set very low values to zero
    r_plot[r_plot <= zero_threshold] <- NA
    # log transform
    r_plot <- log10(r_plot)
    # breaks and legend
    mx <- ceiling(100 * cellStats(r_plot, max)) / 100
    mn <- floor(100 * cellStats(r_plot, min)) / 100
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- sort(c(-2:2, mn, mx))
    lbls <- round(10^lbl_brks, 2)
  } else {
    title <- "Wood Thrush Abundance Uncertainty (SE)"
    # breaks and legend
    mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
    mn <- floor(1000 * cellStats(r_plot, min)) / 1000
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- seq(mn, mx, length.out = 5)
    lbls <- round(lbl_brks, 2)
  }
  
  # abundance
  plot(r_plot, 
       col = pal, breaks = brks, 
       maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE)
  
  # borders
  plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()
  
  # legend
  par(new = TRUE, mar = c(0, 0, 0, 0))
  image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
             smallplot = c(0.25, 0.75, 0.06, 0.09),
             horizontal = TRUE,
             axis.args = list(at = lbl_brks, 
                              labels = lbls,
                              fg = "black", col.axis = "black",
                              cex.axis = 0.75, lwd.ticks = 0.5,
                              padj = -1.5),
             legend.args = list(text = title,
                                side = 3, col = "black",
                                cex = 1, line = 0))
}