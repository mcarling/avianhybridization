ne_land_trim_hz <- st_crop(ne_land, c(ymin=-560000, ymax=1200000, xmin=-5500000, xmax=120000)) # for 8
plot(ne_land_trim_hz)

# transforming to map projection
ne_land_trim_hz <- ne_land_trim_hz %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

# plotting North America
par(mar = c(3.5, 0.25, 0.25, 0.25))
par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(ne_land_trim_final, col = NA, border = NA)
plot(ne_land_trim_final, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# plotting parental taxa abundance ratio
plot(r_plot, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot),
     legend = FALSE, add = TRUE)

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
           #smallplot = c(0.25, 0.75, 0.06, 0.09),
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


# adding lakes and rivers
plot(ne_rivers, col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes, col = "white", border = "#888888", lwd = 0.25, add = TRUE)

# plotting elevation

r_plot_elev <- indigo_r_pred_proj[[3]]

pal_elev <- colorRampPalette(c("#dddddd", "black"))
pal_elev <- pal_elev(26)

# abundance vs. se
#if (nm == "abd") {
#title <- "Ratio P. amoena to P. cyanea"
#title <- expression(paste("Ratio ",italic("P. amoena"), " to ", italic("P. cyanea")))
# set very low values to zero
#r_pred_proj[r_pred_proj <= 0] <- NA ###########
# log transform
#r_plot <- log10(r_plot)

# breaks and legend
mx_elev <- ceiling(100 * cellStats(r_plot_elev, max)) / 100
mn_elev <- floor(100 * cellStats(r_plot_elev, min)) / 100

brks_elev <- seq(mn_elev, mx_elev, length.out = length(pal) + 1)
lbl_brks_elev <- sort(c(-2:2, mn_elev, mx_elev))
lbls_elev <- round(brks_elev, 2)


mx_elev <- ceiling(1000 * cellStats(r_plot_elev, max)) / 1000
mn_elev <- floor(1000 * cellStats(r_plot_elev, min)) / 1000
brks_elev <- seq(mn_elev, mx_elev, length.out = length(pal_elev) + 1)
lbl_brks_elev <- seq(mn_elev, mx_elev, length.out = 5)
lbls_elev <- round(lbl_brks_elev, 2)

# plotting elevation
plot(r_plot_elev, 
     col = pal_elev, breaks = brks_elev, 
     maxpixels = ncell(r_plot_elev),
     legend = FALSE, add = TRUE)

ne_land <- read_sf("data_eBird_best_practices/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()


# plotting
plot(ne_land_trim_final, col = "#dddddd", border = "#888888", lwd = 0.5)
plot(ne_rivers_final, col = "white", lwd = 0.5, add = TRUE)
plot(ne_lakes_final, col = "white", border = "#888888", lwd = 0.25, add = TRUE)





