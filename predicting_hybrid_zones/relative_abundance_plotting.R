## This script uses the predicted relative abundance spatial objects generated
## in the "relative_abundance_modeling.R" script to plot relative 
## abundance of parental taxa and hybrids across North America. It then 
## uses the predicted relative abundance distributions of the two parental species
## to identify areas where hybridization should be possible. Finally, it estimates 
## the relative abundances for each parental species in the hybrid zone.

## As with the other scripts in this folder, much of this code has been taken from 
## the "Best Practices for Using eBird Data" guide by Strimas-Mackey et al. (2020):

## Strimas-Mackey, M., W.M. Hochachka, V. Ruiz-Gutierrez, O.J. Robinson, E.T. Miller,
## T. Auer, S. Kelling, D. Fink, A. Johnston. 2020. Best Practices for Using eBird 
## Data. Version 1.0. https://cornelllabofornithology.github.io/ebird-best-practices/. 
## Cornell Lab of Ornithology, Ithaca, New York. https://doi.org/10.5281/zenodo.3620739


## loading required packages and setting working directory -----------------
library(raster)
library(sf)
library(viridis)
library(fields)
library(tidyverse)

setwd("/Volumes/commons/CarlingLab/eBird Data/Data for looking at relative abundance")


# preparing dataframes with predicted relative abundance values for plotting ------------------

# loading predicted data
indigo_pred <- read.csv("indigo_pred_final.csv", header = TRUE) # for Indigos
lazuli_pred <- read.csv("lazuli_pred_final.csv", header = TRUE) # for Lazulis
hybrids_pred <- read.csv("hybrids_pred_final.csv", header = TRUE) # for hybrids

### converting object to spatial dataframe
# for Indigos
indigo_r_pred <- indigo_pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se, elevation_median) %>% # saving elevation data for plotting below
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
indigo_r_pred <- indigo_r_pred[[-1]]

# for Lazulis
lazuli_r_pred <- lazuli_pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
lazuli_r_pred <- lazuli_r_pred[[-1]]

# for hybrid
hybrid_r_pred <- hybrid_pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
hybrid_r_pred <- hybrid_r_pred[[-1]]


## converting spatial object to same projection as base maps
map_proj <- st_crs(102003)

indigo_r_pred_proj <- projectRaster(indigo_r_pred, crs = map_proj$proj4string, method = "ngb") # for Indigos
lazuli_r_pred_proj <- projectRaster(lazuli_r_pred, crs = map_proj$proj4string, method = "ngb") # for Lazulis
hybrid_r_pred_proj <- projectRaster(hybrid_r_pred, crs = map_proj$proj4string, method = "ngb") # for hybrid

# trimming spatial objects to same limits as base map
e <- extent(-2500000, 2150000, -930000, 1500000) #ymin, ymax, xmin, xmax
indigo_r_pred_proj <- crop(indigo_r_pred_proj, e)
lazuli_r_pred_proj <- crop(lazuli_r_pred_proj, e)
hybrid_r_pred_proj <- crop(hybrid_r_pred_proj, e)


# preparing gis data necessary for making maps ------------------

## load gis data for making maps
map_proj <- st_crs(102003)

ne_land <- read_sf("spatial_data/gis_data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_rivers <- read_sf("spatial_data/gis_data.gpkg", "ne_rivers") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_lakes <- read_sf("spatial_data/gis_data.gpkg", "ne_lakes") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

# trimming spatial objects to same limits as base map
ne_land <- st_crop(ne_land, c(ymin=-930000, ymax=1500000, xmin=-2500000, xmax=2150000))
ne_rivers <- st_crop(ne_rivers, c(ymin=-930000, ymax=1500000, xmin=-2500000, xmax=2150000))
ne_lakes <- st_crop(ne_lakes, c(ymin=-930000, ymax=1500000, xmin=-2500000, xmax=2150000))

## preparing elevation data
r_plot_elev <- indigo_r_pred_proj[["elevation_median"]]

# defining color palette and breaks for shading elevation on maps
pal_elev <- colorRampPalette(c("#dddddd", "black"))
pal_elev <- pal_elev(26)

mx_elev <- ceiling(1000 * cellStats(r_plot_elev, max)) / 1000
mn_elev <- floor(1000 * cellStats(r_plot_elev, min)) / 1000
brks_elev <- seq(mn_elev, mx_elev, length.out = length(pal_elev) + 1)
lbl_brks_elev <- seq(mn_elev, mx_elev, length.out = 5)
lbls_elev <- round(lbl_brks_elev, 2)


# making the relative abundance maps for each species ------------------

### plotting relative abundance for Indigo Bunting
zero_threshold_parentals <- 0.05
zero_threshold_hybrids <- 0.00075 # hybrids are typically much rarer, 
  # need a lower threshold for any predictions to show up on the map

r_plot_indigo <- indigo_r_pred_proj[["abd"]]

par(mar = c(3.5, 0.25, 0.25, 0.25))
# plotting base map
plot(ne_land, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# adding shading for elevation
plot(r_plot_elev, 
     col = pal_elev, breaks = brks_elev, 
     maxpixels = ncell(r_plot_elev),
     legend = FALSE, add = TRUE)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

title <- expression(paste(italic("Passerina cyanea"), " relative abundance"))

# set very low values to zero
r_plot_indigo[r_plot_indigo <= zero_threshold_parentals] <- NA

# log transform
r_plot_indigo <- log10(r_plot_indigo)

# breaks and legend
mx <- ceiling(100 * cellStats(r_plot_indigo, max)) / 100
mn <- floor(100 * cellStats(r_plot_indigo, min)) / 100
brks <- seq(mn, mx, length.out = length(pal) + 1)
lbl_brks <- sort(c(-2:2, mn, mx))
lbls <- round(10^lbl_brks, 2)

# plotting relative abundance
plot(r_plot_indigo, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot_indigo),
     legend = FALSE, add = TRUE)

# adding lakes and rivers
plot(ne_rivers, maxpixels = ncell(r_plot), col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes, maxpixels = ncell(r_plot), col = "white", border = "#888888", lwd = 0.25, add = TRUE)

# adding the legend
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



### plotting relative abundance for Lazuli Bunting
r_plot_lazuli <- lazuli_r_pred_proj[["abd"]]

par(mar = c(3.5, 0.25, 0.25, 0.25))
# plotting base map
plot(ne_land, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# adding shading for elevation
plot(r_plot_elev, 
     col = pal_elev, breaks = brks_elev, 
     maxpixels = ncell(r_plot_elev),
     legend = FALSE, add = TRUE)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

title <- expression(paste(italic("Passerina amoena"), " relative abundance"))

# set very low values to zero
r_plot_lazuli[r_plot_lazuli <= zero_threshold_parentals] <- NA

# log transform
r_plot_lazuli <- log10(r_plot_lazuli)

# breaks and legend
mx <- ceiling(100 * cellStats(r_plot_lazuli, max)) / 100
mn <- floor(100 * cellStats(r_plot_lazuli, min)) / 100
brks <- seq(mn, mx, length.out = length(pal) + 1)
lbl_brks <- sort(c(-2:2, mn, mx))
lbls <- round(10^lbl_brks, 2)

# plotting relative abundance
plot(r_plot_lazuli, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot_lazuli),
     legend = FALSE, add = TRUE)

# adding lakes and rivers
plot(ne_rivers, maxpixels = ncell(r_plot), col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes, maxpixels = ncell(r_plot), col = "white", border = "#888888", lwd = 0.25, add = TRUE)

# adding the legend
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


### plotting relative abundance for hybrid buntings
r_plot_hybrid <- hybrid_r_pred_proj[["abd"]]

par(mar = c(3.5, 0.25, 0.25, 0.25))
# plotting base map
plot(ne_land, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# adding shading for elevation
plot(r_plot_elev, 
     col = pal_elev, breaks = brks_elev, 
     maxpixels = ncell(r_plot_elev),
     legend = FALSE, add = TRUE)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

title <- expression(paste(italic("Passerina amoena"), " x ", italic("cyanea"), " relative abundance"))

# set very low values to zero
r_plot_hybrid[r_plot_hybrid <= zero_threshold_hybrids] <- NA

# log transform
r_plot_hybrid <- log10(r_plot_hybrid)

# breaks and legend
mx <- ceiling(100 * cellStats(r_plot_hybrid, max)) / 100
mn <- floor(100 * cellStats(r_plot_hybrid, min)) / 100
brks <- seq(mn, mx, length.out = length(pal) + 1)
lbl_brks <- sort(c(-2:2, mn, mx))
lbls <- round(10^lbl_brks, 2)

# plotting relative abundance
plot(r_plot_hybrid, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot_hybrid),
     legend = FALSE, add = TRUE)

# adding lakes and rivers
plot(ne_rivers, maxpixels = ncell(r_plot), col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes, maxpixels = ncell(r_plot), col = "white", border = "#888888", lwd = 0.25, add = TRUE)

# adding the legend
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
