## This script uses the filtered and zero-filled eBird datasets prepared in the 
## filtering script, and eBird checklist habitat and bioblim datasets prepared 
## in the hybrid prediction script to run models that determine how different
## variables influence checklist observation count for parental taxa (and also 
## of hybrids if sufficient data exists).
## It then uses the habitat and bioclim prediction surfaces to predict the relative 
## abudnances of parental taxa across North America. Then, it uses these predicted 
## distributions to identify areas where hybridization is possible.

## As with the other scripts in this project, much of this code has been taken from 
## the "Best Practices for Using eBird Data" guide by Strimas-Mackey et al. (2020):

## Strimas-Mackey, M., W.M. Hochachka, V. Ruiz-Gutierrez, O.J. Robinson, E.T. Miller,
## T. Auer, S. Kelling, D. Fink, A. Johnston. 2020. Best Practices for Using eBird 
## Data. Version 1.0. https://cornelllabofornithology.github.io/ebird-best-practices/. 
## Cornell Lab of Ornithology, Ithaca, New York. https://doi.org/10.5281/zenodo.3620739


## loading required packages and setting working directory -----------------
library(raster)
library(mgcv)
library(sf)
library(viridis)
library(fields)
library(dggridR)
library(lubridate)
library(tidyverse)

setwd("/Volumes/commons/CarlingLab/eBird Data/Data for looking at relative abundance")

# set random number seed to ensure fully repeatable results
set.seed(1)


## loading prediction surface and required shape file ------------------------------------------
# reading in the prediction surface
pred_surface <- read_csv("pland-elev-bioclim_prediction-surface_final.csv")
pred_surface <- na.omit(pred_surface)

# latest year of landcover data
max_lc_year <- pred_surface$year[1]

# loading a raster of the prediction surface with bioclimatic and habitat variables at each site
r <- raster("prediction-surface_final.tif") 

# generate hexagonal grid with ~ 5 km betweeen cells for spatiotemporal subsampling
dggs <- dgconstruct(spacing = 5)


## loading and preparing records for each species --------------------------
ebird_data <- c("indigo_bunting_pred_mx_half_july.csv", 
                "lazuli_bunting_pred_mx_half_july.csv", 
                "hybrid_bunting_pred_mx_half_july.csv")

# MODIS habitat covariates, elevation, and bioclim data for all checklist locations
habitat <- read_csv("ebird_pland-elev-bioclim_location-year_final.csv") %>%
  mutate(year = as.integer(year))



ebird_gam <- list()
for(i in ebird_data){
  
  ## reading in the data and set-up
  
  # reading in eBird data
  ebird <- read_csv(i) %>%
    mutate(protocol_type = factor(protocol_type, levels = c("Stationary" , "Traveling"))) %>%
    # remove observations with no count
    filter(!is.na(observation_count)) 
  # combine ebird and habitat data
  ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))
  # need to get rid of NA's
  ebird_habitat <- na.omit(ebird_habitat)
  
  
  ### spatiotemporal subsampling (to mitigate bias in the eBird dataset)
  # determining hexagonal cell id and week number for each checklist
  checklist_cell <- ebird_habitat %>% 
    mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
           week = week(observation_date))
  
  if(sum(checklist_cell$species_observed)>100){
    # sampling one detection and one non-detection checklist per grid cell per week
    ebird_ss <- checklist_cell %>% 
      group_by(species_observed, year, week, cell) %>% 
      sample_n(size = 1) %>% 
      ungroup() %>% 
      select(-cell, -week)}
  else{
    # NOTE: if you try to run a model for hybrids, for which there are likely extremely
    # few observtions, consider retaining all checklists with detections:
    # first divide checklists into whether or not the species was observed, then
    # subsample only for non-detections
    ebird_habitat_detection <- ebird_habitat %>%
      filter(species_observed == TRUE)
    
    ebird_habitat_non_detection <- ebird_habitat %>%
      filter(species_observed == FALSE)
    checklist_cell_nd <- ebird_habitat_non_detection %>%
      mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
             week = week(observation_date))
    
    ebird_ss_nd <- checklist_cell_nd %>%
      group_by(year, week, cell) %>%
      sample_n(size = 1) %>%
      ungroup() %>%
      select(-cell, -week)
    
    ebird_ss <- rbind(ebird_ss_nd, ebird_habitat_detection)
  }
  
  
  ### splitting the subsampled dataset into a training dataset and a 
  ### testing dataset (for evaluating the models)
  
  # Choose which habitat covariates and bioclimatic variables to include in the 
  # model. Consider what habitat types focal taxa prefer and actively avoid.
  hab.bioclim_covs <- c(
    "pland_00_water", 
    "pland_01_evergreen_needleleaf", 
    #"pland_02_evergreen_broadleaf", 
    #"pland_03_deciduous_needleleaf", 
    "pland_04_deciduous_broadleaf", 
    #"pland_05_mixed_forest",
    "pland_06_closed_shrubland", 
    "pland_07_open_shrubland", 
    "pland_08_woody_savanna",#)
    #"pland_09_savanna", 
    "pland_10_grassland", 
    #"pland_11_wetland", 
    "pland_12_cropland", 
    "pland_13_urban", 
    #"pland_14_mosiac", 
    #"pland_15_barren",
    "bio10_mean_temp_warmest_quarter", 
    "bio12_annual_precip",
    "bio18_precip_warmest_quarter"
  )
  
  ebird_split <- ebird_ss %>% 
    # select only the columns to be used in the model
    select(observation_count,
           # effort covariates
           day_of_year, time_observations_started, duration_minutes,
           effort_distance_km, number_observers, protocol_type,
           # habitat and bioclim covariates
           hab.bioclim_covs,
           # also adding latitude, longitude, and elevation
           latitude, longitude, elevation_median)
  
  # split 80/20
  ebird_split <- ebird_split %>% 
    split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
  map_int(ebird_split, nrow)
  
  ebird_gam[[i]] <- ebird_split
}


# running gams to model abundance for each species ------------------------

## gam parameters
# degrees of freedom for smoothing
k <- 5
# degrees of freedom for cyclic time of day smooth
k_time <- 7

# continuous predictors (excluding time, which will be treated as cyclic)
continuous_covs <- ebird_gam[[1]]$train %>% 
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

# trying to do the model with different error distributions
# explicitly specify where the knots should occur for time_observations_started
# this ensures that the cyclic spline joins the variable at midnight
# this won't happen by default if there are no data near midnight
time_knots <- list(time_observations_started = seq(0, 24, length.out = k_time))

# negative binomial
m_nb <- gam(gam_formula,
            data = ebird_gam[[i]]$train, 
            family = "nb",
            knots = time_knots)



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

indigo_m_nb <- m_nb
pred_model <- m_nb

seq_tod <- seq(0, 24, length.out = 300)
tod_df <- ebird_split$train %>% 
  # find average pland habitat covariates
  select(starts_with("pland"), starts_with("bio"), latitude, longitude, elevation_median) %>% 
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
       title = "Effect of observation start time on Indigo Bunting reporting",
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
#### maybe need to unscale latitude, longitude, and elevation before fitting model to prediction surface!!!
##### if you do try to scale latitude, longitude, and elevation, maybe try to scale those variables 
##### in the prediction dataframe, then predict, and then unscale in the new dataframe for mapping
indigo_pred <- predict(pred_model, newdata = pred_surface_eff, # for Indigos
#lazuli_pred <- predict(pred_model, newdata = pred_surface_eff, # for Lazulis
                       #hybrid_pred <- predict(pred_model, newdata = pred_surface_eff, # for hybrids                       
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

write.csv(indigo_pred, "indigo_pred_nb_final_final.csv", na = "", row.names=FALSE) # for Indigos


## for model validation, look at distribution of residuals

## also looking at influential observations








# making species distribution maps ----------------------------------------
# loading predicted data
indigo_pred <- read.csv("indigo_pred_nb_final.csv", header = TRUE) # for Indigos
lazuli_pred <- read.csv("lazuli_pred_nb_final.csv", header = TRUE) 

# for Indigos
indigo_r_pred <- indigo_pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
indigo_r_pred <- indigo_r_pred[[-1]]

# load gis data for making maps
map_proj <- st_crs(102003)
ne_land <- read_sf("data_eBird_best_practices/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()


#filtering down to new datasets
ne_land_trim_final <- st_crop(ne_land, c(ymin=-930000, ymax=1500000, xmin=-2500000, xmax=2150000))

# now converting to native modis projection
ne_land_modis <- ne_land_trim_final %>% 
  # project to the native modis projection
  st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))
plot(ne_land_modis)




ne_land_trim2 <- ne_land_trim %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf("data_eBird_best_practices/gis-data.gpkg", "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data_eBird_best_practices/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf("data_eBird_best_practices/gis-data.gpkg", "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()







# plotting distirbutions for the two parental species ---------------------
zero_threshold <- 0.05

#for (nm in names(indigo_r_pred)) {
r_plot <- r_pred_proj_lazuli[["abd"]]

par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(ne_land, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

# abundance vs. se
#if (nm == "abd") {
title <- "Predicted relative abundance"
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
# } else {
#   title <- "hybrid bunting Abundance Uncertainty (SE)"
#   # breaks and legend
#   mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
#   mn <- floor(1000 * cellStats(r_plot, min)) / 1000
#   brks <- seq(mn, mx, length.out = length(pal) + 1)
#   lbl_brks <- seq(mn, mx, length.out = 5)
#   lbls <- round(lbl_brks, 2)
# }

# abundance
plot(r_plot, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot),
     legend = FALSE, add = TRUE)

# borders
#plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
#box()

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



# plotting hybrid zone ----------------------------------------------------
indigo_pred_abd <- indigo_pred %>%
  select(latitude, longitude, abd) %>%
  rename(indigo_abd = abd)
lazuli_pred_abd <- lazuli_pred %>%
  select(latitude, longitude, abd) %>%
  rename(lazuli_abd = abd)

indigo_pred_abd$lazuli_abd <- lazuli_pred_abd$lazuli_abd
hybrid_pred_abd <- indigo_pred_abd
#hybrid_pred_abd <- inner_join(indigo_pred_abd, lazuli_pred_abd, by = c("latitude", "longitude"))
#hybrid_pred_abd <- merge(indigo_pred_abd, lazuli_pred_abd)

hybrid_pred_abd <- hybrid_pred_abd %>%
  #mutate(hybrid_zone = ifelse(indigo_abd >= 0.1 & lazuli_abd >= 0.1, TRUE, FALSE))
  mutate(hybrid_zone = ifelse(indigo_abd >= 0.05 & lazuli_abd >= 0.05, TRUE, FALSE))
#mutate(hybrid_zone = ifelse(indigo_abd >= 0.025 & lazuli_abd >= 0.025, TRUE, FALSE))
#mutate(hybrid_zone = ifelse(indigo_abd >= 0.01 & lazuli_abd >= 0.01, TRUE, FALSE))

# looking at the relative proprotions of yes's and no's
length(hybrid_pred_abd$hybrid_zone)
sum(hybrid_pred_abd$hybrid_zone)


### converting object to spatial dataframe
hybrid_r_pred_abd <- hybrid_pred_abd %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(hybrid_zone) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
hybrid_r_pred_abd <- hybrid_r_pred_abd[[-1]]


hybrid_pred_abd_overlap <- hybrid_pred_abd %>%
  filter(hybrid_zone == TRUE)

dim(lazuli_pred_abd)
dim(indigo_pred_abd)
dim(hybrid_pred_abd)
setdiff(indigo_pred_abd$latitude, lazuli_pred_abd$latitude) 
test <- indigo_pred_abd$latitude %in% lazuli_pred_abd$latitude
unique(test)
length(test)
sum(test)
lazuli_pred_abd %>%
  
  25.19792, 25.19792
names(indigo_pred)

# project predictions
map_proj <- st_crs(102003)
r_pred_proj_hz <- projectRaster(hybrid_r_pred_abd, crs = map_proj$proj4string, method = "ngb")



par(mfrow = c(1, 1))
#for (nm in names(indigo_r_pred)) {
r_plot <- r_pred_proj_hz[[1]]

par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(ne_land, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

# abundance vs. se
#if (nm == "abd") {
title <- "P. cyanea x P. amoena hybrid zone"
# set very low values to zero
r_pred_proj_hz[r_pred_proj_hz <= 0] <- NA
# log transform
#r_plot <- log10(r_plot)
# breaks and legend
mx <- 1
mn <- 0
brks <- seq(mn, mx, length.out = length(pal) + 1)
lbl_brks <- sort(c(-2:2, mn, mx))
lbls <- round(10^lbl_brks, 2)
# } else {
#   title <- "Indigo Bunting Abundance Uncertainty (SE)"
#   # breaks and legend
#   mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
#   mn <- floor(1000 * cellStats(r_plot, min)) / 1000
#   brks <- seq(mn, mx, length.out = length(pal) + 1)
#   lbl_brks <- seq(mn, mx, length.out = 5)
#   lbls <- round(lbl_brks, 2)
# }

# abundance
plot(r_pred_proj_hz[[1]], 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot),
     legend = FALSE, add = TRUE)

# borders
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


# Plotting hybrid zone shaded by proporitonal abundance of parenta --------
# looking at proportion of each parental species in each place
hybrid_pred_abd <- hybrid_pred_abd %>%
  mutate(lazuli_prop = (lazuli_abd / (indigo_abd + lazuli_abd))) %>%
  mutate(lazuli_diff = (lazuli_abd - indigo_abd))

# subsetting hybrid data to include only data in area of overlap
hybrid_pred_abd_overlap <- hybrid_pred_abd %>%
  filter(hybrid_zone == TRUE)


### converting object to spatial dataframe
hybrid_r_pred_abd_overlap <- hybrid_pred_abd_overlap %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(hybrid_zone, lazuli_prop, lazuli_diff) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
hybrid_r_pred_abd_overlap <- hybrid_r_pred_abd_overlap[[-1]]


## plotting area of overlap and proportion of Lazulis
map_proj <- st_crs(102003)
r_pred_proj <- projectRaster(hybrid_r_pred_abd_overlap, crs = map_proj$proj4string, method = "ngb")

par(mfrow = c(1, 1))
#for (nm in names(indigo_r_pred)) {
r_plot <- r_pred_proj[["lazuli_prop"]]

par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(ne_land, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# modified plasma palette
# plasma_rev <- rev(plasma(25, end = 0.9))
# gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
# pal <- c(gray_int(4)[2], plasma_rev)
# length(pal)
pal <- colorRampPalette(c("red", "blue"))
pal <- pal(26)

# abundance vs. se
#if (nm == "abd") {
title <- "Proportion P. amoena"
# set very low values to zero
r_pred_proj[r_pred_proj <= 0] <- NA ###########
# log transform
#r_plot <- log10(r_plot)
# breaks and legend
mx <- 1
mn <- 0
brks <- seq(mn, mx, length.out = length(pal) + 1)
lbl_brks <- sort(c(-2:2, mn, mx))
lbls <- round(brks, 1)
lbls2 <- lbls[seq(1, length(lbls), 2)]
brks2 <- seq(0,1, length = 11)
lbls2 <- seq(0,1, length = 11)
# } else {
#   title <- "Indigo Bunting Abundance Uncertainty (SE)"
#   # breaks and legend
#   mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
#   mn <- floor(1000 * cellStats(r_plot, min)) / 1000
#   brks <- seq(mn, mx, length.out = length(pal) + 1)
#   lbl_brks <- seq(mn, mx, length.out = 5)
#   lbls <- round(lbl_brks, 2)
# }


# abundance
plot(r_pred_proj[["lazuli_prop"]], 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot),
     legend = FALSE, add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
           smallplot = c(0.25, 0.75, 0.06, 0.09),
           horizontal = TRUE,
           axis.args = list(at = brks2, 
                            labels = lbls2,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))
}
