#loading required packages
library(maps)
library(maptools)
library(sf)
library(tmap)
library(tmaptools)
library(mapproj)
library(auk)
library(tidyverse)
library(graphics)

#=============================================================================================
#     mapping hybridization rates by county based on the raw eBird data
#=============================================================================================

#filtering the entire dataset down to May and June and to RBSA, RNSA, and hybrids
ebd_dir <- "/pfs/tsfs1/gscratch/pdoughe1" #if doing on Teton, just make sure 
  #you use working directory where eBird txt file is stored
ebd_dir <- "/Volumes/commons/CarlingLab/eBird Data/" #if doing this on personal 
  #machine connected to PetaLibrary

f_in_ebd <- file.path(ebd_dir, "ebd_US_relJan-2020.txt")
ebd_filters <- auk_ebd(f_in_ebd) %>% 
  auk_species(species = c("Red-naped Sapsucker", "Red-breasted Sapsucker", "Red-naped x Red-breasted Sapsucker (hybrid)")) %>%
  auk_state(state = c("US-CA", "US-OR", "US-WA")) %>% #maybe also include British Columbia, Alberta, 
  #Idaho, Colorado, Arizona, New Mexico, Nevada, Texas, Mexico?
  auk_date(date = c("*-05-09", "*-06-05")) #the date range over which that Shawn collected in 2012
#this currently looks at all years

f_out_ebd <- "sapsuckers.txt"
if (!file.exists(f_out_ebd)) { 
  ebd_filtered <- auk_filter(ebd_filters, file = f_out_ebd, overwrite = TRUE)
} #the if statement prevents this code from running if the filtered file already exists

#reading in the data into R
if (!exists("sapsuckers")) { 
  sapsuckers <- read_ebd("sapsuckers.txt", rollup=FALSE)
}


###preparing the eBird dataset for analyses

#getting rid of X's
sapsuckers <- sapsuckers %>% 
  filter(observation_count !="X")
#converting observation count to a numeric vector
sapsuckers$observation_count = as.numeric(sapsuckers$observation_count)


###calculating proportion of hybrids in each county

#making a vector of counties that appear in the dataset
sapsucker_counties <- unique(sapsuckers$county_code)

#calculating the rate of hybridization in each county (# of (rows of) hybrids/ # of (rows of) all birds)
# county_proportions <- sapply(sapsucker_counties, function(w,x,y){
#   x <- sapsuckers %>%
#     filter(county_code == w)
#   y <- x %>%
#     filter(category=="hybrid")
#   nrow(y)/nrow(x)
# })

#calculating the rate of hybridization in each county (# of hybrids/ # of all birds)
county_proportions <- sapply(sapsucker_counties, function(w,x,y){
  x <- sapsuckers %>%
    filter(county_code == w)
  y <- x %>% 
    filter(category=="hybrid")
  (sum(y$observation_count))/(sum(x$observation_count))
})
typeof(county_proportions)
county_proportions_df <- as.data.frame(county_proportions)

#renaming the columns of the dataframe
names(county_proportions_df)[1] <- "hybrid_proportion"
county_proportions_df <- tibble::rownames_to_column(county_proportions_df, var = "county_code")


###mapping county hybrid proportions

#making a dataframe of unique counties with their full name and state
states_counties <- unique(sapsuckers[c("state", "county", "county_code")])

#binding the two dataframes
states_counties <- na.omit(states_counties)
county_proportions_df <- na.omit(county_proportions_df)
county_prop_states <- inner_join(states_counties, 
                                 county_proportions_df, 
                                 by="county_code")

#loading the shape file
load(url("https://github.com/mgimond/ES218/blob/gh-pages/Data/counties48.RData?raw=true"))
head(cnty)

#creating ID column in eBird data file to match county shape file
county_prop_states <- county_prop_states %>%
  mutate(ID = paste(tolower(state), tolower(county), sep = ","))

#finding out which county names are in the hybrid proportion dataset and not the shapefile
setdiff(county_prop_states$ID, cnty$ID)

#joining the county geometry shape file and hybrid proportion data
county_prop_states_map <- inner_join(cnty, county_prop_states, by="ID" )

#creating the map
tm_shape(county_prop_states_map) +
  tm_polygons(col='hybrid_proportion',palette='PuRd',n=7,style='jenks',
  #for style: 'jenks' minimizes variation within groups, 'fisher' maximizes variation among groups
              title='proportion hybrids',border.lwd=0, border.alpha=1) +
  tm_legend(position=c("right","top"),frame=FALSE,outside=TRUE,
            main.title='May 9th - June 5th')


#=============================================================================================
#     mapping hybridization rates by county based on Shawn's data
#=============================================================================================

###preparing the data

setwd("~/Documents/eBird_hybrid_analyses")
#reading in the dataframe
sapsucker_shawn <- read.csv("Sapsucker_localities,h_indices.csv", sep=",", header=T)

#defining what individuals are hybrids and which are species
sapsucker_shawn$category <- ifelse(sapsucker_shawn$Hindex_1 <= 0.1 | 
                                     sapsucker_shawn$Hindex_1 >= 0.9, 
                                   "species", "hybrid")

#converting state and county to character vectors
sapsucker_shawn$state <- as.character(sapsucker_shawn$state)
sapsucker_shawn$county <- as.character(sapsucker_shawn$county)

###calculating the rate of hybridiaztion (# of hybrids / # of all birds) in each county

#making a vector of counties that appear in the dataset
sapsucker_shawn_counties <- unique(sapsucker_shawn$county)

##defining a function that calculates the rate of hybridization in each county (# of hybrids/ # of all birds)
county_proportions_shawn <- sapply(sapsucker_shawn_counties, function(w,x,y){
  x <- sapsucker_shawn %>%
    filter(county == w)
  y <- x %>% 
    filter(category=="hybrid")
  nrow(y)/nrow(x)
})
county_proportions_shawn <- sapply(Sapsucker_shawn_counties, county_fun_shawn)
typeof(county_proportions_shawn)
county_proportions_shawn_df <- as.data.frame(county_proportions_shawn)

#renaming the columns of the dataframe
names(county_proportions_shawn_df)[1] <- "hybrid_proportion"
county_proportions_shawn_df <- tibble::rownames_to_column(county_proportions_shawn_df, var = "county")

###mapping the hybridization rate in each county

#making a dataframe of unique counties with their full name and state
states_counties_shawn <- unique(sapsucker_shawn[c("state", "county")])

#binding the two dataframes
states_counties_shawn <- na.omit(states_counties_shawn)
county_proportions_shawn_df <- na.omit(county_proportions_shawn_df)
county_prop_states_shawn <- inner_join(states_counties_shawn, 
                                       county_proportions_shawn_df, 
                                       by="county")

#creating an ID column for each county to match the shape file with geometry of each county
county_prop_states_shawn <- county_prop_states_shawn %>%
  mutate(ID = paste(tolower(state), tolower(county), sep = ","))

#finding out which county names are in the hybrid proportion dataset and not the shapefile
setdiff(county_prop_states_shawn$ID, cnty$ID)

#binding the hybrid rate data and the shapefile
county_prop_states_shawn_map <- inner_join(cnty, county_prop_states_shawn, by="ID" )

# #for creating states shape file (need to download this file)
# shp <- mutate(shp, STFIPS = stringr::str_sub(GEOID, 1, 2))
# states <- shp %>%
#   aggregate_map(by = "STFIPS")
# tm_shape(states) +
#   tm_borders(col = "black")

#creating mappable state borders by aggregating county geometry data
cnty$state <- gsub( ",.*", "", cnty$ID) 
sapsucker_hz <- c("washington", "oregon", "california", "idaho", "nevada") #for whatever reason it doesn't work if you try aggregate all state borders at once
states <- cnty[cnty$state %in% sapsucker_hz,] %>%
  aggregate_map(by = "state")
# tm_shape(states) +
#   tm_borders(col = "black")

#mapping Shawn's hybrid rate data with state borders
tm_shape(county_prop_states_shawn_map) +
  tm_polygons(col='hybrid_proportion',palette='PuRd',n=7,style='fisher',
              title='proportion hybrids',border.lwd=0, border.alpha=1) +
  tm_shape(states) +
  tm_borders(col = "black") +
  tm_legend(position=c("right","top"),frame=FALSE,outside=TRUE,
            main.title='Shawn, 2012')


#=============================================================================================
# mapping hybridization rates by county based on the eBird data and Shawn's data at the same time
#=============================================================================================

###merging Shawn's data with eBird in a single dataframe
#filtering down eBird data to just counties that Shawn sampled from
county_prop_states_map_reduced <- county_prop_states_map %>%
  filter(geometry %in% county_prop_states_shawn_map$geometry)

##mapping the eBird data in the counties that Shawn sampled from
# tm_shape(county_prop_states_map_reduced) +
#   tm_polygons(col='hybrid_proportion',palette='PuRd',n=7,style='fisher',
#               title='hybrid proportion',border.lwd=0, border.alpha=1) +
#   tm_shape(states) +
#   tm_borders(col = "black") +
#   tm_legend(position=c("right","top"),frame=TRUE,outside=TRUE,
#             main.title='eBird')


###combining the two in a single map

#adding an identifier column to each dataset
county_prop_states_shawn_map$dataset <- "Shawn"
county_prop_states_map_reduced$dataset <- "eBird"
county_prop_states_map_reduced <- county_prop_states_map_reduced %>%
  select(-county_code)

#combining Shawn's data and the eBird data into a single dataframe
all <- rbind(county_prop_states_shawn_map,county_prop_states_map_reduced)

#creating maps with the same legend
tm_shape(all) +
  tm_polygons(col='hybrid_proportion',palette='PuRd',n=9,style='fisher',
              title='Hybrid Proportion',border.lwd=0, border.alpha=1) +
  tm_facets(by = "dataset", free.coords = TRUE, ncol = 2) +
  tm_shape(states) +
  tm_borders(col = "black") +
  tm_legend(position=c("right","top"),frame=FALSE,outside=TRUE,
            main.title='')


###plotting difference between Shawn's data and eBird data in each county

#filtering down each dataframe (I guess they're technically spatial objects now) to just relevant columns
county_prop_states_shawn_diff <- county_prop_states_shawn_map %>%
  select(ID, hybrid_proportion)
county_prop_states_diff <- county_prop_states_map_reduced %>%
  select(ID, hybrid_proportion)

#combining the two dataframes horizontally
county_prop_diff_all <- cbind(county_prop_states_diff,
                                county_prop_states_shawn_diff)

#calculating the difference in hybridization rates between Shawn's data and the eBird data for each county
county_prop_diff_all$difference <- county_prop_diff_all$hybrid_proportion.1 - county_prop_diff_all$hybrid_proportion

#plotting differences in hybridization rates between Shawn's data and the eBird data for each county

ggplot(county_prop_diff_all, aes(difference, ID)) +
  geom_point() +
  geom_vline(xintercept = 0) +
  labs(y="", x = "Difference between Shawn and eBird")


#=============================================================================================
# comparing hybrid proportion calculations between a maximum likelihood and Bayesian framework
#=============================================================================================


#=============================================================================================
#     investigating how variation in eBird effort influences hybrid proportion map
#=============================================================================================

###as a first pass, looking at correlation between hybrid proportion, total checklists submitted, 
  #and between total sapsuckers reported

##defining the hybrid zone as counties where both parental species cooccur or where 
  ##either parental species occurs with hybrids
county_taxa <- list()
for(i in unique(sapsuckers$county_code)){
  county_taxa[[i]] <- unique(subset(sapsuckers,county_code == i)$common_name)
}
##filtering down to counties where at least two taxa have been reported
counties<-lapply(county_taxa, function(x)(length(x)>1))
#a dataframe of counties and whether or not there were two buntings reported in each county
counties_df <- do.call(rbind,counties)
counties_df <- as.data.frame(counties_df)
#renaming columns in this dataframe
names(counties_df)[1] <- "hybrid_chance"
counties_df <- tibble::rownames_to_column(counties_df, var = "county_code")
#making a vector of counties where at least two buntings have been reported
hybrid_counties_df <- counties_df %>%
  filter(hybrid_chance ==TRUE)
#the list of counties where at least two buntings have been recorded
hybrid_counties <- hybrid_counties_df$county_code
#filtering initial dataset to just include counties where at least two sapsuckers have been recorded
hybrid_zone_dataset <- sapsuckers %>%
  filter(county_code %in% hybrid_counties)


# ##now calculating the total number of checklists submitted in each county in the hybrid zone 
# #!!!! NOTE: this is only the number of checklists with sapsuckers!!!!!
# 
# #making a vector of counties that appear in the dataset
# hz_counties <- unique(hybrid_zone_dataset$county_code)
# 
# #performing a fucntion that calculates the total number of checklists submitted in each county of the hybrid zone
# hz_county_checklists <- sapply(hz_counties, function(x,y){
#   y <- hybrid_zone_dataset %>%
#     filter(county_code == x)
#   (length(unique(y$checklist_id)))
# })
# 
# hz_county_checklists_df <- as.data.frame(hz_county_checklists)
# #renaming the columns of the dataframe
# names(hz_county_checklists_df)[1] <- "checklists_submitted"
# hz_county_checklists_df <- tibble::rownames_to_column(hz_county_checklists_df, var = "county_code")


##this calculates the total number of checklists submitted in each county, 
  ##regardless of if they include sapsuckers or not
#filtering the entire dataset down to May and June and to RBSA, RNSA, and hybrids
ebd_dir <- "/pfs/tsfs1/gscratch/pdoughe1" #if doing on Teton, just make sure 
#you use working directory where eBird txt file is stored
ebd_dir <- "/Volumes/commons/CarlingLab/eBird Data/" #if doing this on personal machine

f_in_ebd <- file.path(ebd_dir, "ebd_US_relJan-2020.txt")
ebd_filters <- auk_ebd(f_in_ebd) %>% 
  #auk_species(species = c("Red-naped Sapsucker", "Red-breasted Sapsucker", "Red-naped x Red-breasted Sapsucker (hybrid)")) %>%
  auk_state(state = c("US-CA", "US-OR", "US-WA")) %>% #maybe also include British Columbia, Alberta, 
  #Idaho, Colorado, Arizona, New Mexico, Nevada, Texas, Mexico?
  auk_date(date = c("*-05-09", "*-06-05")) #the date range over which that Shawn collected in 2012
#this currently looks at all years

f_out_ebd <- "WA.OR.CA.txt"
if (!file.exists(f_out_ebd)) { 
  ebd_filtered <- auk_filter(ebd_filters, file = f_out_ebd, overwrite = TRUE)
} #the if statement prevents this code from running if the filtered file already exists

#reading in the data into R
if (!exists("WA.OR.CA")) { 
  WA.OR.CA <- read_ebd("WA.OR.CA.txt", rollup=FALSE)
}

#filtering down to counties where sapsuckers have opportunity to hybridize
hz_full_dataset <- WA.OR.CA %>%
  filter(county_code %in% hybrid_counties)

#calculating the number of checklists submitted to eBird in these counties
#performing a fucntion that calculates the total number of checklists submitted in each county of the hybrid zone
hz_county_checklists <- sapply(hz_counties, function(x,y){
  y <- hz_full_dataset %>%
    filter(county_code == x)
  (length(unique(y$checklist_id)))
})

hz_county_checklists_df <- as.data.frame(hz_county_checklists)
#renaming the columns of the dataframe
names(hz_county_checklists_df)[1] <- "checklists_submitted"
hz_county_checklists_df <- tibble::rownames_to_column(hz_county_checklists_df, var = "county_code")


##calculating total number of sapsuckers reported in each county of the hybrid zone

hz_county_sapsuckers_df<- aggregate(observation_count ~ county_code, data=hybrid_zone_dataset, sum, na.rm=T)
#renaming the columns of the dataframe
names(hz_county_sapsuckers_df)[2] <- "sapsuckers_reported"


##calculating the proportion of sapsuckers reported as hybrids in each county of the hybrid zone
hz_county_props <- sapply(hz_counties, function(w,x,y){
  x <- hybrid_zone_dataset %>%
    filter(county_code == w)
  y <- x %>% 
    filter(category=="hybrid")
  sum(y$observation_count)/sum(x$observation_count)
})

hz_county_props_df <- as.data.frame(hz_county_props)
#renaming the columns of the dataframe
names(hz_county_props)[1] <- "hybrid_proportion"
hz_county_props_df <- tibble::rownames_to_column(hz_county_props_df, var = "county_code")


##now binding all three variables in a dataframe and examining correlation among them
hz_data <- merge(hz_county_checklists_df, hz_county_sapsuckers_df) 
hz_data <- merge(hz_data, hz_county_props_df)

MyVar <- c("checklists_submitted", "sapsuckers_reported", "hz_county_props")
pairs(hz_data[,MyVar], lower.panel=panel.cor)



###creating a variogram for the eBird data

#filtering eBird atatset down to counties where Shawn sampled from (maybe just contiguous counties?)





                       