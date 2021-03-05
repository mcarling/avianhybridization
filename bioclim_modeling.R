## In this script, we model bioclimatic suitability for hybrid 
## Passerina amoena x cyanea based first on museum records, and 
## then based on museum records and eBird records together. We 
## use the dismo package to calculate how similar areas are based 
## on how similar they are to locations where hybrids were reported, 
## based on the 19 bioclimatic variables. Much of this code was inspired 
## by this guide put together by Jeff Oliver: https://jcoliver.github.io/learn-r/011-species-distribution-models.html.

library(sp)
library(raster)
library(maptools)
library(rgdal)
library(dismo)
library(tidyverse)
library(lubridate)
library(sf)

setwd("/Volumes/commons/CarlingLab/eBird Data/Data for looking at relative abundance")


#=============================================================================================
#     preparing bioclim data
#=============================================================================================

# organizing workspace
dir.create(path = "data")
dir.create(path = "output")

# accessing worldclim data
bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5,
                        path = "data/")

## can run models with only bioclim variables we're interested in. Otherwise, models will run based on all 19 variables.
relevant_vars <- c(10,12,18)
bioclim.relevant <- bioclim.data[[relevant_vars]]


#=============================================================================================
#     loading and preparing records of hybrids
#=============================================================================================

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


#=============================================================================================
#     modeling hybrid bioclimatic distribution based on museum records
#=============================================================================================

# Determine geographic extent of our data
max.lat <- ceiling(max(museum.data$decimallatitude))
min.lat <- floor(min(museum.data$decimallatitude))
max.lon <- ceiling(max(museum.data$decimallongitude))
min.lon <- floor(min(museum.data$decimallongitude))
geographic.extent <- extent(x = c(min.lon - 5, max.lon + 5, min.lat - 5, max.lat + 5))

# Crop bioclim data to geographic extent of museum records for hybrids
bioclim.data2 <- crop(x = bioclim.data, y = geographic.extent)
# bioclim.relevant2 <- crop(x = bioclim.relevant, y = geographic.extent)


# Build species distribution model
bc.model <- bioclim(x = bioclim.data2, p = museum.data)
# bc.model <- bioclim(x = bioclim.relevant2, p = museum.data)


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


#=============================================================================================
#     modeling hybrid bioclimatic distribution based on museum records and eBird data
#=============================================================================================

# Determine geographic extent of our data
max.lat <- ceiling(max(museum.ebird.data$decimallatitude))
min.lat <- floor(min(museum.ebird.data$decimallatitude))
max.lon <- ceiling(max(museum.ebird.data$decimallongitude))
min.lon <- floor(min(museum.ebird.data$decimallongitude))
geographic.extent <- extent(x = c(min.lon - 5, max.lon + 5, min.lat - 5, max.lat + 5))

# Crop bioclim data to geographic extent of museum records for hybrids
bioclim.data2 <- crop(x = bioclim.data, y = geographic.extent)
# bioclim.relevant2 <- crop(x = bioclim.relevant, y = geographic.extent)


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
