## In this script, we model bioclimatic suitability for hybrid 
## Passerina amoena x cyanea based first on museum records, and 
## then based on museum records and eBird records together. We 
## use the dismo package to calculate how similar areas are based 
## on how similar they are to locations where hybrids were reported, 
## based on the 19 bioclimatic variables. Much of this code was inspired 
## by this guide put together by Jeff Oliver: https://jcoliver.github.io/learn-r/011-species-distribution-models.html.


# set up ------------------------------------------------------------------

# organizing workspace
dir.create(path = "data")
dir.create(path = "output")

library(sp)
library(raster)
library(maptools)
library(rgdal)
library(dismo)
library(tidyverse)
library(lubridate)
library(sf)

# accessing worldclim data
bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5,
                        path = "data/")

## filtering down bioclim dataset to only the variables we're interested in
relevant_vars <- c(10,12,18)

bioclim.relevant <- bioclim.data[[relevant_vars]]

!!!!!!! also in the modeling script, look at covariance between predictor variables



# laoding and preparing data --------------------------------------------

#### reading in and preparing museum data
museum.data <- read.csv("Passerina_VertNet_hybrids.csv")
# Check the data to make sure it loaded correctly
summary(museum.data)

# filtering down to just lat and long columns, reordering them for model
museum.data <- museum.data[, c("decimallongitude", "decimallatitude")]

# drop entries that don't have lat long available
museum.data <- museum.data[!is.na(museum.data$decimallatitude), ]


#### reading in and preparing eBird records for hybrids
hybrid_bunting_pred <- read.csv("hybrid_bunting_pred_mx_half_july.csv", header= TRUE)
# preparing eBird data
hybrid_present_ebird <- hybrid_bunting_pred %>%
  filter(species_observed == TRUE) %>%
  mutate(decimallongitude = longitude,
         decimallatitude = latitude)

# filtering down to just lat and long columns, reordering them for model
ebird.data <- hybrid_present_ebird[, c("decimallongitude", "decimallatitude")]


#### combining ebird and museum data
museum.ebird.data <- rbind(museum.data, ebird.data)




# # Load the data to use for our base map
# data(wrld_simpl)
# data(state)
# #plot(state)
# 
# # Plot the base map
# plot(wrld_simpl, 
#      xlim = c(min.lon - 5, max.lon + 5),
#      ylim = c(min.lat - 5, max.lat + 5),
#      axes = TRUE, 
#      col = "grey95")


#plot(wrld_simpl, add = TRUE)

# Add the points for individual observation
# points(x = museum.ebird.data$decimallongitude, 
#        y = museum.ebird.data$decimallatitude, 
#        col = "green", 
#        pch = 20, 
#        cex = 0.75)
# points(x = museum.data$decimallongitude, 
#        y = museum.data$decimallatitude, 
#        col = "red", 
#        pch = 20, 
#        cex = 0.75)
# # And draw a little box around the graph
# box()


# modeling bioclimatic distribution based on museum records ---------------

# Determine geographic extent of our data
max.lat <- ceiling(max(museum.data$decimallatitude))
min.lat <- floor(min(museum.data$decimallatitude))
max.lon <- ceiling(max(museum.data$decimallongitude))
min.lon <- floor(min(museum.data$decimallongitude))
geographic.extent <- extent(x = c(min.lon - 5, max.lon + 5, min.lat - 5, max.lat + 5))

# Crop bioclim data to geographic extent of museum records for hybrids
bioclim.data2 <- crop(x = bioclim.data, y = geographic.extent)
bioclim.relevant2 <- crop(x = bioclim.relevant, y = geographic.extent)


# Build species distribution model
bc.model <- bioclim(x = bioclim.data2, p = museum.data)
bc.model <- bioclim(x = bioclim.relevant2, p = museum.data)


# Predict presence from model
predict.presence <- dismo::predict(object = bc.model, 
                                   x = bioclim.data2, 
                                   ext = geographic.extent)
# need to specify the predict function in the dismo package, 
# instead of the predict function in the raster package


## now plotting to have the same background map and color
## as the other figures

# projecting to match the map of north america
predict.presence2 <- projectRaster(predict.presence, crs = map_proj$proj4string, method = "ngb") # for hybrid

# cropping to same geogrpahic area as other maps
e <- extent(-2500000, 2150000, -930000, 1500000) #ymin, ymax, xmin, xmac
predict.presence2 <- crop(predict.presence2, e)

# any areas with predicted suitabilty below this threshold are set to zero
zero_threshold <- 0.05 # the same number as the abundance threshold that we used for hybrids earlier

par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(ne_land_trim_final, col = NA, border = NA)
plot(ne_land_trim_final, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# adding a layer for elevation (from another script)
plot(r_plot_elev, 
     col = pal_elev, breaks = brks_elev, 
     maxpixels = ncell(r_plot_elev),
     legend = FALSE, add = TRUE)


# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

# abundance
title <- expression(paste("Potential bioclimatic distribution of ", italic("Passerina amoena"), " x ", italic("cyanea")))
# set very low values to zero
predict.presence2[predict.presence2 <= zero_threshold] <- NA

# breaks and legend
mx <- ceiling(1000 * cellStats(predict.presence2, max)) / 1000
mn <- floor(1000 * cellStats(predict.presence2, min)) / 1000
brks <- seq(0, 1, length.out = length(pal) + 1)
lbl_brks <- seq(0, 1, length.out = 5)
lbls <- round(lbl_brks, 2)


# abundance
plot(predict.presence2, 
     col = pal, breaks = brks, 
     maxpixels = ncell(predict.presence2),
     legend = FALSE, add = TRUE)

# # Add original observations
# points(museum.data$decimallongitude, museum.data$decimallatitude, col = "olivedrab", pch = 20, cex = 0.75)
# box()

plot(ne_rivers, maxpixels = ncell(r_plot), col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes, maxpixels = ncell(r_plot), col = "white", border = "#888888", lwd = 0.25, add = TRUE)

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



# modeling bioclimatic distribution based on museum records and eBird data ---------------

# Determine geographic extent of our data
max.lat <- ceiling(max(museum.ebird.data$decimallatitude))
min.lat <- floor(min(museum.ebird.data$decimallatitude))
max.lon <- ceiling(max(museum.ebird.data$decimallongitude))
min.lon <- floor(min(museum.ebird.data$decimallongitude))
geographic.extent <- extent(x = c(min.lon - 5, max.lon + 5, min.lat - 5, max.lat + 5))

# Crop bioclim data to geographic extent of museum records for hybrids
bioclim.data2 <- crop(x = bioclim.data, y = geographic.extent)
bioclim.relevant2 <- crop(x = bioclim.relevant, y = geographic.extent)


# Build species distribution model
bc.model <- bioclim(x = bioclim.data2, p = museum.ebird.data)

# Predict presence from model
predict.presence <- dismo::predict(object = bc.model, 
                                   x = bioclim.data2, 
                                   ext = geographic.extent)

## plotting
# projecting to match the map of north america
predict.presence3 <- projectRaster(predict.presence, crs = map_proj$proj4string, method = "ngb") # for hybrid

# cropping to same geogrpahic area
e <- extent(-2500000, 2150000, -930000, 1500000) #ymin, ymax, xmin, xmac
predict.presence3 <- crop(predict.presence3, e)

# setting zero threshold
zero_threshold <- 0.00075

par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(ne_land_trim_final, col = NA, border = NA)
plot(ne_land_trim_final, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# adding elevation
plot(r_plot_elev, 
     col = pal_elev, breaks = brks_elev, 
     maxpixels = ncell(r_plot_elev),
     legend = FALSE, add = TRUE)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

# abundance
title <- expression(paste("Potential bioclimatic distribution of ", italic("Passerina amoena"), " x ", italic("cyanea")))
# set very low values to zero
predict.presence3[predict.presence3 <= zero_threshold] <- NA

# breaks and legend
mx <- ceiling(1000 * cellStats(predict.presence3, max)) / 1000
mn <- floor(1000 * cellStats(predict.presence3, min)) / 1000
brks <- seq(0, 1, length.out = length(pal) + 1)
lbl_brks <- seq(0, 1, length.out = 5)
lbls <- round(lbl_brks, 2)


# abundance
plot(predict.presence3, 
     col = pal, breaks = brks, 
     maxpixels = ncell(predict.presence3),
     legend = FALSE, add = TRUE)

# # Add original observations
# points(museum.data$decimallongitude, museum.data$decimallatitude, col = "olivedrab", pch = 20, cex = 0.75)
# box()

plot(ne_rivers, maxpixels = ncell(r_plot), col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes, maxpixels = ncell(r_plot), col = "white", border = "#888888", lwd = 0.25, add = TRUE)

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














# Use the bioclim data files for sampling resolution
bil.files <- list.files(path = "data/wc2-5", 
                        pattern = "*.bil$", 
                        full.names = TRUE)

# We only need one file, so use the first one in the list of .bil files
mask <- raster(bil.files[1])

#Randomly sample points (same number as our observed points)
background <- randomPoints(mask = mask,     # Provides resolution of sampling points
                           n = nrow(museum.data),      # Number of random points
                           ext = geographic.extent, # Spatially restricts sampling
                           extf = 1.25)             # Expands sampling a little bit

head(background)


# Plot the base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95",
     main = "Presence and pseudo-absence points")

# Add the background points
points(background, col = "grey30", pch = 1, cex = 0.75)

# Add the observations
points(x = museum.data$decimallongitude, 
       y = museum.data$decimallatitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)

box()


# Arbitrarily assign group 1 as the testing data group
testing.group <- 1

# Create vector of group memberships
group.presence <- kfold(x = museum.data, k = 5) # kfold is in dismo package


head(group.presence)

# Should see even representation in each group
table(group.presence)


# Separate observations into training and testing groups
presence.train <- museum.data[group.presence != testing.group, ]
presence.test <- museum.data[group.presence == testing.group, ]

# Repeat the process for pseudo-absence points
group.background <- kfold(x = background, k = 5)
background.train <- background[group.background != testing.group, ]
background.test <- background[group.background == testing.group, ]


# Build a model using training data
bc.model <- bioclim(x = bioclim.data, p = presence.train)

# Predict presence from model (same as previously, but with the update model)
predict.presence <- dismo::predict(object = bc.model, 
                                   x = bioclim.data, 
                                   ext = geographic.extent)


# Use testing data for model evaluation
bc.eval <- evaluate(p = presence.test,   # The presence testing data
                    a = background.test, # The absence testing data
                    model = bc.model,    # The model we are evaluating
                    x = bioclim.data)    # Climatic variables for use by model

# Determine minimum threshold for "presence"
bc.threshold <- threshold(x = bc.eval, stat = "spec_sens")


# Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

# Only plot areas where probability of occurrence is greater than the threshold
plot(predict.presence > bc.threshold, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "olivedrab"))

# And add those observations
points(x = museum.data$decimallongitude, 
       y = museum.data$decimallatitude, 
       col = "black",
       pch = "+", 
       cex = 0.75)

# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()











##### combining Bioclim with MODIS and eBird data ------------------------------------------
# extracting bioclim values for eBird checklists --------------------------
#converting observation date to date, not a factor
lazuli_bunting_pred <- read.csv("lazuli_bunting_pred_mx_half_july.csv", header = TRUE) # for Lazuli
lazuli_bunting_pred$observation_date <- as_date(lazuli_bunting_pred$observation_date)

indigo_bunting_pred <- read.csv("indigo_bunting_pred_mx_half_july.csv", header = TRUE) # for Indigo
indigo_bunting_pred$observation_date <- as_date(indigo_bunting_pred$observation_date)

hybrid_bunting_pred <- read.csv("hybrid_bunting_pred_mx_half_july.csv", header = TRUE) # for hybrid
hybrid_bunting_pred$observation_date <- as_date(hybrid_bunting_pred$observation_date)


#lazuli_sf <- lazuli_bunting_pred %>% # for Lazuli
#indigo_sf <- indigo_bunting_pred %>% # for Indigo
hybrid_sf <- hybrid_bunting_pred %>% # for hybrids
  # filtering down to unique localities, don't care about year for this, I don't think
  # we also don't need to create a buffer around each point, all we need is the value for each point
  distinct(locality_id, latitude, longitude) %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% # have to change crs? no, you shouldn't need to
  # transform to bioclim projection
  st_transform(crs = projection(bioclim.data))

# need to pull out each variable from the bioclim dataset to add to lazuli sf dataset!!!!!!!
#lazuli_locs <- st_set_geometry(lazuli_sf, NULL) %>% # for Lazulis
#indigo_locs <- st_set_geometry(indigo_sf, NULL) %>% # for Lazulis
hybrid_locs <- st_set_geometry(hybrid_sf, NULL) %>% # for Lazulis
  mutate(id = row_number())

relevant_vars <- c(10,12,18)

#bioclim_values = vector('list', nlayers(bioclim.data))
bioclim_values = vector('list', length(relevant_vars))

for (i in relevant_vars){
  # for Lazuli
  #bioclim_checklists <- raster::extract(bioclim.data[[i]], lazuli_sf)
  #bioclim_values[[i]] <- tibble(locality_id = lazuli_locs$locality_id, bioclim_var = names(bioclim.data[[i]]), bioclim_value = bioclim_checklists)

  # for Indigo
  #bioclim_checklists <- raster::extract(bioclim.data[[i]], indigo_sf)
  #bioclim_values[[i]] <- tibble(locality_id = indigo_locs$locality_id, bioclim_var = names(bioclim.data[[i]]), bioclim_value = bioclim_checklists)
  
  # for hybrid
  bioclim_checklists <- raster::extract(bioclim.data[[i]], hybrid_sf)
  bioclim_values[[i]] <- tibble(locality_id = hybrid_locs$locality_id, bioclim_var = names(bioclim.data[[i]]), bioclim_value = bioclim_checklists)
  
}

bioclim_values_all <- do.call(rbind, bioclim_values)

# renaming bioclim variable names

# creating a dataframe of bioclim variable names
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


# combining with new names
#lazuli_bioclim <- bioclim_values_all %>% # for Lazuli
#indigo_bioclim <- bioclim_values_all %>% # for indigo
hybrid_bioclim <- bioclim_values_all %>% # for hybrid
  inner_join(bioclim_names, by = "bioclim_var") %>% 
  arrange(bioclim_var) %>% 
  select(-bioclim_var)

#lazuli_bioclim <- lazuli_bioclim %>% # for Lazuli
#indigo_bioclim <- indigo_bioclim %>% # for Indigo
hybrid_bioclim <- hybrid_bioclim %>% # for hybrid
  pivot_wider(names_from = bioclim_name, 
              values_from = bioclim_value, 
              values_fill = list(bioclim_value = NA))

# saving dataframes, but probably don't even need to save different 
# dataframes for different species, because they should be the same, right?
# for Lazuli
write_csv(lazuli_bioclim, "lazuli_bioclim_final.csv", na = "")
lazuli_bioclim <- read.csv("lazuli_bioclim_final.csv")

# for Indigo
write_csv(indigo_bioclim, "indigo_bioclim_final.csv", na = "")
indigo_bioclim <- read.csv("indigo_bioclim_final.csv")

# for hybrid
write_csv(hybrid_bioclim, "hybrid_bioclim_final.csv", na = "")
hybrid_bioclim <- read.csv("hybrid_bioclim_final.csv")


# extracting bioclim values for prediction surface ------------------------

# adding bioclim values to prediction surface

# including progress bar
library(progress)
pb <- progress_bar$new(
  format = "  downloading [:bar] :percent eta: :eta",
  total = 100, clear = FALSE, width= 60)


relevant_vars <- c(10,12,18) # only looking at relevant bioclim variables
# bioclim_values_pred = vector('list', nlayers(bioclim.data))
bioclim_values_pred = vector('list', length(relevant_vars))

for (i in relevant_vars){ # maybe should try rcells instead just to be consistent?
  bioclim_centers <- raster::extract(bioclim.data[[i]], r_centers)  # check to make sure projections are the same!!!
  bioclim_values_pred[[i]] <- tibble(id = r_centers$id, bioclim_var = names(bioclim.data[[i]]), bioclim_value = bioclim_centers)
  #bioclim_values[[i]] <- tibble(locs$locality_id, names(bioclim.data[[i]]), bioclim_checklists)
  # r_centers2 <- r_centers
  # r_centers2$bioclim <- bioclim_values
  # bioclim_values_pred[[i]] <- r_centers2
  # bioclim_values_pred[[i]]$bioclim_var <- names(bioclim.data[[i]])
  pb$tick(100/3) # keeping track of progress
}

bioclim_values_pred_all <- do.call(rbind, bioclim_values_pred)

unique(bioclim_values_pred_all$bioclim_var) # making sure it worked

# combining with new names and spreading
bioclim_values_pred_all <- bioclim_values_pred_all %>%
  inner_join(bioclim_names, by = "bioclim_var") %>% 
  arrange(bioclim_var) %>% 
  select(-bioclim_var)

head(bioclim_values_pred_all)

#bioclim_values_pred_all_test <- spread(bioclim_values_pred_all, bioclim_name, bioclim)

bioclim_values_pred_all <- bioclim_values_pred_all %>%
  pivot_wider(names_from = bioclim_name, 
              values_from = bioclim_value, 
              values_fill = list(bioclim_value = NA))


# join in coordinates
bioclim_pred <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(bioclim_values_pred_all, by = "id")

write_csv(bioclim_pred, "pred_bioclim_r_centers.csv", na = "")
