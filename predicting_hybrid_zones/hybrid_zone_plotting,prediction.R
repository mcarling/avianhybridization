## This script uses the predicted relative abundance spatial objects generated
## in the "relative_abundance_modeling.R" script to first estimate where hybridization 
## occurs between two focal taxa, and then estimates the abundance rations of these 
## two parental taxa throughout the hybrid zone.

## loading required packages and setting working directory -----------------
library(raster)
library(sf)
library(viridis)
library(fields)
library(tidyverse)

setwd("/Volumes/commons/CarlingLab/eBird Data/Data for looking at relative abundance")


# defining the hybrid zone and preparing a spatial object for mapping ------------------------------------

# prediction surface with estimated relative abundance for each parental species
indigo_pred <- read.csv("indigo_pred_nb_final.csv")
lazuli_pred <- read.csv("lazuli_pred_nb_final.csv")

# merging the two prediction surfaces
indigo_pred_abd <- indigo_pred %>%
  select(latitude, longitude, abd) %>%
  rename(indigo_abd = abd)
lazuli_pred_abd <- lazuli_pred %>%
  #select(latitude, longitude, abd) %>%
  select(abd) %>%
  rename(lazuli_abd = abd)

#hybrid_pred_abd <- merge(indigo_pred_abd, lazuli_pred_abd) # this probably better, if you can get it to work
hybrid_pred_abd <- cbind(indigo_pred_abd, lazuli_pred_abd)

## plotting areas in which hybridization is possible, that being where both
## parental taxa occur above a threshold abundance level. Choose the abundance threshold that you think is 
## most appropriate. You can adjust the threshold after looking at maps. You can also have different thresholds 
## for each parental species.
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

### converting hybrid zone object to spatial dataframe
hybrid_r_pred_abd <- hybrid_pred_abd_overlap %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  #select(hybrid_zone, lazuli_rat, indigo_abd) %>% 
  select(hybrid_zone) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
hybrid_r_pred_abd <- hybrid_r_pred_abd[[-1]]


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
r_plot_elev <- indigo_r_pred_proj[[3]]

# defining color palette and breaks for shading elevation on maps
pal_elev <- colorRampPalette(c("#dddddd", "black"))
pal_elev <- pal_elev(26)

mx_elev <- ceiling(1000 * cellStats(r_plot_elev, max)) / 1000
mn_elev <- floor(1000 * cellStats(r_plot_elev, min)) / 1000
brks_elev <- seq(mn_elev, mx_elev, length.out = length(pal_elev) + 1)
lbl_brks_elev <- seq(mn_elev, mx_elev, length.out = 5)
lbls_elev <- round(lbl_brks_elev, 2)


# plotting the hybrid zone ----------------------------------------------------

### plotting ratio of Lazulis to Indigos in areas of overlap
map_proj <- st_crs(102003)
hz_pred_proj <- projectRaster(hybrid_r_pred_abd, crs = map_proj$proj4string, method = "ngb")

r_plot_hz <- hz_pred_proj[["lazuli_rat"]]

par(mar = c(3.5, 0.25, 0.25, 0.25))
# plotting base map
plot(ne_land, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# adding shading for elevation
plot(r_plot_elev, 
     col = pal_elev, breaks = brks_elev, 
     maxpixels = ncell(r_plot_elev),
     legend = FALSE, add = TRUE)

# color palette for shading the hybrid zone
pal <- colorRampPalette(c("red", "blue")) # make this a solid colr (e.g., purple, if you just 
 # want to plot the hybrid zone, not shaded by the ratio of parentals
pal <- pal(26)

title <- expression(paste("Ratio ",italic("P. amoena"), " to ", italic("P. cyanea")))

# log transform
r_plot_hz <- log10(r_plot_hz)

# breaks and legend
mx <- ceiling(1000 * cellStats(r_plot_hz, max)) / 1000
mn <- floor(1000 * cellStats(r_plot_hz, min)) / 1000
brks <- seq(mn, mx, length.out = length(pal) + 1)
lbl_brks <- seq(mn, mx, length.out = 5)
lbls <- round(lbl_brks, 2)

 
# plotting hybrid zone shaded by abundance ratio of the two parental taxa
plot(r_plot_hz, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot_hz),
     legend = FALSE, add = TRUE)

# adding bodies of water for reference
plot(ne_rivers, col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes, col = "white", border = "#888888", lwd = 0.25, add = TRUE)
#plot(ne_land_trim, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

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



# plotting predicted relative abundance of each parental taxon at a specific latitude --------------------

### for the Niobrara River, which goes through the hybrid zone

# filtering down to approximate latitude of the Niobrara
hybrid_pred_abd_niobrara <- hybrid_pred_abd_overlap %>%
  #filter(latitude >= 41.77662 & latitude <= 42.78812)
  filter(latitude >= 42.78 & latitude <= 42.78812) # depending on how wide your 
   # latitudinal band is, may need to calculate the mean relative abudnance for 
   # each longitudinal point using the agregate function
#filter(latitude == 42.78125)
unique(hybrid_pred_abd_niobrara$latitude)

# exploratory plots to see appropriate plot dimensions
plot(lazuli_abd ~ longitude, data = hybrid_pred_abd_niobrara, type = "l", 
     xlab = "longitude", ylab = "relative abundance", col = "red")
plot(indigo_abd ~ longitude, data = hybrid_pred_abd_niobrara, type = "l", 
     xlab = "longitude", ylab = "relative abundance", col = "blue")


## making a new plot from scratch
quartz.options(height=8,width=10,dp=72) # making a new plotting window, specifying height, width, resolution
plot.new();
par(mfcol= c(1,1));
frame();
plot.window(xlim=c(-103.5, -100.5), ylim=c(0,0.40)); # setting dimensions of plotting window

lines(lazuli_abd ~ longitude, data = hybrid_pred_abd_niobrara, col = "red", lwd = 2)
lines(lazuli_lcl ~ longitude, data = hybrid_pred_abd_niobrara, col = "red", lty = "dashed")
lines(lazuli_pcl ~ longitude, data = hybrid_pred_abd_niobrara, col = "red", lty = "dashed")
lines(indigo_abd ~ longitude, data = hybrid_pred_abd_niobrara, col = "blue", lwd = 2)
lines(indigo_lcl ~ longitude, data = hybrid_pred_abd_niobrara, col = "blue", lty = "dashed")
lines(indigo_pcl ~ longitude, data = hybrid_pred_abd_niobrara, col = "blue", lty = "dashed")

#use mtext() to add labels to axis.
axis(side=1, at=seq(-103.5, -100.5, by = 0.5), cex.axis=1);
axis(side=2, at=seq(0,0.5, by=0.1), cex.axis = 1, las=1); # at argument tells you where to put tickmarks


mtext("Predicted Relative Abundance ", side = 2, cex=1.5, line=2.7,par(family=""))

mtext("Longitude", side = 1, cex=1.5, line=2.5,family="sans",font=1)

legend("topright", legend=c(expression(italic("P. amoena")),
                       expression(italic("P. cyanea"))), 
       col=c("red", "blue"), lwd = 2, cex=1.25)


###for the Platte River, which also goes through the hybrid zone

# filtering down to approximate latitude of the Niobrara
hybrid_pred_abd_platte <- hybrid_pred_abd_overlap %>%
  filter(latitude >= 41.01041 & latitude <= 41.01043)
unique(hybrid_pred_abd_platte$latitude)

# exploratory plots to see appropriate plot dimensions
plot(lazuli_abd ~ longitude, data = hybrid_pred_abd_platte, type = "l", 
     xlab = "longitude", ylab = "relative abundance", col = "red")
plot(indigo_abd ~ longitude, data = hybrid_pred_abd_platte, type = "l", 
     xlab = "longitude", ylab = "relative abundance", col = "blue")

# making a nice plot from scratch
quartz.options(height=8,width=10,dp=72) ##making a new plotting window, specifying height, width, resolution
plot.new();
par(mfcol= c(1,1));
frame();
plot.window(xlim=c(-104, -100.5), ylim=c(0,0.40)); ##setting dimensions of plotting window

lines(lazuli_abd ~ longitude, data = hybrid_pred_abd_platte, col = "red", lwd = 2)
lines(lazuli_lcl ~ longitude, data = hybrid_pred_abd_platte, col = "red", lty = "dashed")
lines(lazuli_pcl ~ longitude, data = hybrid_pred_abd_platte, col = "red", lty = "dashed")
lines(indigo_abd ~ longitude, data = hybrid_pred_abd_platte, col = "blue", lwd = 2)
lines(indigo_lcl ~ longitude, data = hybrid_pred_abd_platte, col = "blue", lty = "dashed")
lines(indigo_pcl ~ longitude, data = hybrid_pred_abd_platte, col = "blue", lty = "dashed")

#use mtext() to add labels to axis.
axis(side=1, at=seq(-104, -100.5, by = 0.5), cex.axis=1);
axis(side=2, at=seq(0,0.5, by=0.1), cex.axis = 1, las=1); # at argument tells you where to put tickmarks


mtext("Predicted Relative Abundance ", side = 2, cex=1.5, line=2.7,par(family=""))

mtext("Longitude", side = 1, cex=1.5, line=2.5,family="sans",font=1)

legend("topright", legend=c(expression(italic("P. amoena")),
                            expression(italic("P. cyanea"))), 
       col=c("red", "blue"), lwd = 2, cex=1.25)
