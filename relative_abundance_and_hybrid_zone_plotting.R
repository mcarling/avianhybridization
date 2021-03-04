## this script uses the relative abundance spatial objects to plot relative 
## abundance maps for the two species as well as areas of hybridization

## loading required packages and setting working directory -----------------
library(raster)
library(mgcv)
library(viridis)
library(fields)
library(dggridR)
library(tidyverse)




# getting rid of predicted abundance in water for indigo ------------------

indigo_pred_full <- merge(indigo_pred, pred_surface)


# converting areas that are totally surrounded by water to 0 predicted abundance for Indigo Buntings
#indigo_pred_full <- transform(indigo_pred_full, abd = ifelse(pland_00_water==1, 0, abd))

indigo_pred <- indigo_pred_full %>% 
  select(latitude, longitude, abd, abd_se, abd_lcl, abd_ucl, elevation_median)

names(indigo_pred)

# converting prediction dataframes to spatial objects for making maps --------

# for Indigos
indigo_r_pred <- indigo_pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se, elevation_median) %>% 
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

# for hybrids
hybrid_r_pred <- hybrid_pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
hybrid_r_pred <- hybrid_r_pred[[-1]]

# converting raster to same projection as map spatial object
map_proj <- st_crs(102003)
indigo_r_pred_proj <- projectRaster(indigo_r_pred, crs = map_proj$proj4string, method = "ngb") # for Indigos
lazuli_r_pred_proj <- projectRaster(lazuli_r_pred, crs = map_proj$proj4string, method = "ngb") # for Indigos
hybrid_r_pred_proj <- projectRaster(hybrid_r_pred, crs = map_proj$proj4string, method = "ngb") # for hybrid

e <- extent(-2500000, 2150000, -930000, 1500000) #ymin, ymax, xmin, xmac
indigo_r_pred_proj <- crop(indigo_r_pred_proj, e)
lazuli_r_pred_proj <- crop(lazuli_r_pred_proj, e)
hybrid_r_pred_proj <- crop(hybrid_r_pred_proj, e)

ne_rivers <- st_crop(ne_rivers, c(ymin=-930000, ymax=1500000, xmin=-2500000, xmax=2150000))
ne_lakes <- st_crop(ne_lakes, c(ymin=-930000, ymax=1500000, xmin=-2500000, xmax=2150000))


# making relative abundance maps for the two parental species -------------

# any expected abundances below this threshold are set to zero
zero_threshold <- 0.05
zero_threshold <- 0.00075 # for hybrids

r_plot_elev <- indigo_r_pred_proj[[3]]
r_plot_indigo <- indigo_r_pred_proj[[1]]
r_plot_lazuli <- lazuli_r_pred_proj[[1]]
r_plot_hybrid <- hybrid_r_pred_proj[[1]]


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
title <- expression(paste(italic("Passerina cyanea"), " relative abundance"))
title <- expression(paste(italic("Passerina amoena"), " relative abundance"))
title <- expression(paste(italic("Passerina amoena"), " x ", italic("cyanea"), " relative abundance"))
# set very low values to zero
r_plot_indigo[r_plot_indigo <= zero_threshold] <- NA
r_plot_lazuli[r_plot_lazuli <= zero_threshold] <- NA
r_plot_hybrid[r_plot_hybrid <= zero_threshold] <- NA

# log transform
r_plot_indigo <- log10(r_plot_indigo)
r_plot_lazuli <- log10(r_plot_lazuli)
r_plot_hybrid <- log10(r_plot_hybrid)
# breaks and legend
mx <- ceiling(100 * cellStats(r_plot_indigo, max)) / 100
mn <- floor(100 * cellStats(r_plot_indigo, min)) / 100
brks <- seq(mn, mx, length.out = length(pal) + 1)
lbl_brks <- sort(c(-2:2, mn, mx))
lbls <- round(10^lbl_brks, 2)


# abundance
plot(r_plot_indigo, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot),
     legend = FALSE, add = TRUE)


plot(ne_rivers, maxpixels = ncell(r_plot), col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes, maxpixels = ncell(r_plot), col = "white", border = "#888888", lwd = 0.25, add = TRUE)



# borders
#plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
#plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
#plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
#box()

########## for rivers and lakes
ne_rivers <- ne_download(scale = 50, type = 'rivers_lake_centerlines', 
                         category = 'physical',
                         returnclass = "sf") %>%
  st_set_precision(1e6) %>%
  st_union()

ne_lakes <- ne_download(scale = 50, type = 'lakes', 
                        category = 'physical',
                        returnclass = "sf") %>%
  st_set_precision(1e6) %>%
  st_union()

ne_rivers <- ne_rivers %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_lakes <- ne_lakes %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

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



# plotting hybrid zone ----------------------------------------------------

# plotting areas in which hybridization is possible, that being where both
# parental taxa occur above a threshold abundance level

#prediction surface with estimated relative abundances of Indigos and Lazuli Buntings
indigo_pred <- read.csv("indigo_pred_nb_final.csv")
lazuli_pred <- read.csv("lazuli_pred_nb_final.csv")

#prediction surface with estimated relative abundances of Indigos and Lazuli Buntings
#prediction surface with estimated relative abundances of Indigos and Lazuli Buntings
indigo_pred_abd <- indigo_pred %>%
  select(latitude, longitude, abd) %>%
  rename(indigo_abd = abd)
lazuli_pred_abd <- lazuli_pred %>%
  #select(latitude, longitude, abd) %>%
  select(abd) %>%
  rename(lazuli_abd = abd)

#hybrid_pred_abd <- merge(indigo_pred_abd, lazuli_pred_abd)
hybrid_pred_abd <- cbind(indigo_pred_abd, lazuli_pred_abd)

hybrid_pred_abd <- hybrid_pred_abd %>%
  #mutate(hybrid_zone = ifelse(indigo_abd >= 0.1 & lazuli_abd >= 0.1, TRUE, FALSE))
  #mutate(hybrid_zone = ifelse(indigo_abd >= 0.075 & lazuli_abd >= 0.075, TRUE, FALSE))
mutate(hybrid_zone = ifelse(indigo_abd >= 0.05 & lazuli_abd >= 0.05, TRUE, FALSE))
# mutate(hybrid_zone = ifelse(indigo_abd >= 0.025 & lazuli_abd >= 0.025, TRUE, FALSE))
  #mutate(hybrid_zone = ifelse(indigo_abd >= 0.01 & lazuli_abd >= 0.025, TRUE, FALSE))
  
# looking at the relative proprotions of cells in hybrid zone vs. outside of hybrid zone
length(hybrid_pred_abd$hybrid_zone)
sum(hybrid_pred_abd$hybrid_zone)

# subsetting hybrid data to include only data in area of overlap
hybrid_pred_abd_overlap <- hybrid_pred_abd %>%
  filter(hybrid_zone == TRUE)

# looking at proportion of each parental species in each place
hybrid_pred_abd_overlap <- hybrid_pred_abd_overlap %>%
  #mutate(lazuli_prop = (lazuli_abd / (indigo_abd + lazuli_abd))) %>%
  #mutate(lazuli_diff = (lazuli_abd - indigo_abd)) %>%
  mutate(lazuli_rat = (lazuli_abd / indigo_abd))


### converting object to spatial dataframe
hybrid_r_pred_abd <- hybrid_pred_abd_overlap %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  #select(hybrid_zone, lazuli_rat, indigo_abd) %>% 
  select(hybrid_zone) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
hybrid_r_pred_abd <- hybrid_r_pred_abd[[-1]]


## plotting ratio of Lazulis to Indigos in areas of overlap
map_proj <- st_crs(102003)
r_pred_proj <- projectRaster(hybrid_r_pred_abd, crs = map_proj$proj4string, method = "ngb")

r_plot <- r_pred_proj[["lazuli_rat"]]
r_plot <- r_pred_proj

# plotting North America
par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(ne_land_trim_final, col = NA, border = NA)
plot(ne_land_trim_final, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# modified plasma palette
# plasma_rev <- rev(plasma(25, end = 0.9))
# gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
# pal <- c(gray_int(4)[2], plasma_rev)
# length(pal)
pal <- colorRampPalette(c("blue", "purple"))
pal <- pal(26)

# abundance vs. se
#if (nm == "abd") {
#title <- "Ratio P. amoena to P. cyanea"
title <- expression(paste("Ratio ",italic("P. amoena"), " to ", italic("P. cyanea")))
# set very low values to zero
#r_pred_proj[r_pred_proj <= 0] <- NA ###########
# log transform
r_plot <- log10(r_plot)

# breaks and legend
mx <- ceiling(100 * cellStats(r_plot, max)) / 100
mn <- floor(100 * cellStats(r_plot, min)) / 100

# brks <- seq(mn, mx, length.out = length(pal) + 1)
# lbl_brks <- sort(c(-2:2, mn, mx))
# lbls <- round(brks, 2)
# lbls2 <- lbls[seq(1, length(lbls), 2)]
# brks2 <- seq(0,1, length = 11)
# lbls2 <- seq(0,1, length = 11)
# } else {
#   title <- "Indigo Bunting Abundance Uncertainty (SE)"
#   # breaks and legend
   mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
   mn <- floor(1000 * cellStats(r_plot, min)) / 1000
   brks <- seq(mn, mx, length.out = length(pal) + 1)
   lbl_brks <- seq(mn, mx, length.out = 5)
   lbls <- round(lbl_brks, 2)
# }

 
# abundance
plot(r_plot, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot),
     legend = FALSE, add = TRUE)

# plot(r_plot, 
#      col = pal, breaks = brks, 
#      maxpixels = ncell(r_plot),
#      legend = FALSE, add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
box()



plot(ne_rivers, col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes, col = "white", border = "#888888", lwd = 0.25, add = TRUE)
plot(ne_land_trim, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)


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











# project predictions
map_proj <- st_crs(102003)
# r_pred_proj <- projectRaster(indigo_r_pred, crs = map_proj$proj4string, method = "ngb") # for Indigos
r_pred_proj <- projectRaster(lazuli_r_pred, crs = map_proj$proj4string, method = "ngb") # for Lazulis
# r_pred_proj <- projectRaster(hybrid_r_pred, crs = map_proj$proj4string, method = "ngb") # for hybrid

r_pred_proj_lazuli

# load gis data for making maps
map_proj <- st_crs(102003)
ne_land <- read_sf("data_eBird_best_practices/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

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



par(mfrow = c(2, 1))
for (nm in names(lazuli_r_pred)) {
  r_plot <- r_pred_proj_lazuli[[nm]]
  
  par(mar = c(3.5, 0.25, 0.25, 0.25))
  # set up plot area
  plot(ne_land, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  
  # modified plasma palette
  plasma_rev <- rev(plasma(25, end = 0.9))
  gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
  pal <- c(gray_int(4)[2], plasma_rev)
  
  # abundance vs. se
  if (nm == "abd") {
    title <- "hybrid bunting relative abundance"
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
    title <- "hybrid bunting abundance uncertainty (SE)"
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





