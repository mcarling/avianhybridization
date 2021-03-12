### In this script, we use the eBird datasets prepared in the "all_hybridization_rates" script to 
### calculate hybridization rates for a few specific hybrid systems, rather than the entire eBird 
### dataset. We demonstrate how imposing different geographic filters can greatly influence estimates 
### of system-specific hybridization rates.

library(auk)
library(tidyverse)

# Set-up and preparing the eBird data ------------------------------------------------------------------

setwd('/Volumes/commons/CarlingLab/eBird Data')
# setwd('/pfs/tsfs1/gscratch/pdoughe1') # if you're running this script on Teton

#reading in the eBird data for each year separately
ebd2010 <- read_ebd("ebd_filtered2010.txt", rollup=FALSE)
ebd2011 <- read_ebd("ebd_filtered2011.txt", rollup=FALSE)
ebd2012 <- read_ebd("ebd_filtered2012.txt", rollup=FALSE)
ebd2013 <- read_ebd("ebd_filtered2013.txt", rollup=FALSE)
ebd2014 <- read_ebd("ebd_filtered2014.txt", rollup=FALSE)
ebd2015 <- read_ebd("ebd_filtered2015.txt", rollup=FALSE)
ebd2016 <- read_ebd("ebd_filtered2016.txt", rollup=FALSE)
ebd2017 <- read_ebd("ebd_filtered2017.txt", rollup=FALSE)
ebd2018 <- read_ebd("ebd_filtered2018.txt", rollup=FALSE)

# Combining the eBird data for each year into a list:
ebd10to18 <- list(ebd2010, ebd2011, ebd2012, ebd2013, ebd2014, ebd2015, ebd2016, ebd2017, ebd2018)


#making a dataframe of the frequency of all hybrid combinations that occur in the dataset
  #this step isn't necessary, just tells you the number of each hybrid combination to help you 
  #decide which ones to include in the analyses
hybrids <- unique(ebd2018[ebd2018$category=="hybrid",]$common_name)
hybrid_count <- function(x,y){
  y <- ebd2018 %>% # can do for other years
    filter(common_name == x)
  (sum((y$observation_count)))
}
hybrid_counts <- sapply(hybrids, hybrid_count)
typeof(hybrid_counts)
hybrid_counts_df <- as.data.frame(hybrid_counts)


#defining a character vector for each hybridizing group that we want to include
buntings <- c("Indigo Bunting", "Lazuli Bunting", 
              "Lazuli x Indigo Bunting (hybrid)")
orioles <- c("Baltimore Oriole", "Bullock's Oriole", 
             "Bullock's x Baltimore Oriole (hybrid)")
sapsuckers <- c("Red-naped Sapsucker", "Red-breasted Sapsucker", 
                "Red-naped x Red-breasted Sapsucker (hybrid)")
chickadees <- c("Black-capped Chickadee", "Carolina Chickadee", 
                "Carolina x Black-capped Chickadee (hybrid)" )
titmice <- c("Tufted Titmouse", "Black-crested Titmouse", 
             "Tufted x Black-crested Titmouse (hybrid)")
v.warblers <- c("Blue-winged Warbler", "Golden-winged Warbler", 
                "Golden-winged x Blue-winged Warbler (hybrid) ", "Brewster's Warbler (hybrid)", "Lawrence's Warbler (hybrid)")
s.warblers <- c("Townsend's Warbler", "Hermit Warbler", 
                "Townsend's x Hermit Warbler (hybrid)")
grosbeaks <- c("Rose-breasted Grosbeak", "Black-headed Grosbeak", 
               "Rose-breasted x Black-headed Grosbeak (hybrid)") #there are very few of these
towhees <- c("Eastern Towhee", "Spotted Towhee", 
             "Spotted x Eastern Towhee (hybrid)") #there are very few of these
gulls <- c("Western Gull", "Glaucous-winged Gull", 
           "Western x Glaucous-winged Gull (hybrid)")
ducks1 <- c("Mallard", "Mottled Duck", 
            "Mallard x Mottled Duck (hybrid)")
ducks2 <- c("Mallard", "American Black Duck",
            "American Black Duck x Mallard (hybrid)")
sparrows <- c("Saltmarsh Sparrow", "Nelson's Sparrow", 
              "Nelson's x Saltmarsh Sparrow (hybrid)")
grebes <- c("Western Grebe", "Clark's Grebe", 
            "Western x Clark's Grebe (hybrid)")
ibises <- c("Glossy Ibis", "White-faced Ibis", 
            "Glossy x White-faced Ibis (hybrid)")
quail <- c("California Quail", "Gambel's Quail", 
           "California x Gambel's Quail (hybrid) ") #only one of these in 2016

#combining each vector into a single list
hybrid_systems <- list(buntings, orioles, sapsuckers, chickadees, 
                       titmice, v.warblers, s.warblers, 
                       grosbeaks, towhees, gulls, ducks1, ducks2, 
                       sparrows, grebes, ibises, quail)


#calculating 3 different hybridization rates for hybrid systems of choice for each year in the ebird dataset
for(year in ebd10to18)
{
  year <- year %>% 
    filter(observation_count !="X") # this gets rid of X's
  year$observation_count = as.numeric(year$observation_count)
  
  #printing the year
  year$ebdyear <- (sub("-.*", "", year$observation_date))
  print(unique(year$ebdyear))
  
  for(hybrid_system in hybrid_systems) #running these steps for each defined hybrid system
  {filtered_dataset <- year %>%
    filter(common_name %in% hybrid_system)
  
  ###1) hybridization rates across the entire range of the two focal species
  hybrids_rangewide <- filtered_dataset %>% 
    filter(category=="hybrid")
  parental1_rangewide <-filtered_dataset %>% 
    filter(common_name == hybrid_system[1])
  parental2_rangewide <- filtered_dataset %>% 
    filter(common_name == hybrid_system[2])
  total_rangewide <- filtered_dataset %>% 
    filter(category=="hybrid"| category=="species")
  
  total_rate_rangewide <- (sum(hybrids_rangewide$observation_count))/(sum(total_rangewide$observation_count))
  parental1_rate_rangewide <- (sum(hybrids_rangewide$observation_count))/(sum(parental1_rangewide$observation_count))
  parental2_rate_rangewide <- (sum(hybrids_rangewide$observation_count))/(sum(parental2_rangewide$observation_count))
  
  
  ###2)hybridization rates in counties where both parental species occur, or where either 
  #parental species occurs with hybrids
  #defining a for loop that generates a list of unique species in each county
  county_taxa <- list()
  for(i in unique(filtered_dataset$county_code)){
    county_taxa[[i]] <- unique(subset(filtered_dataset,county_code == i)$common_name)
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
  #filtering initial dataset tojust include counties where at least two buntings have been recorded
  hybrid_zone_dataset <- filtered_dataset %>%
    filter(county_code %in% hybrid_counties)
  
  #calculating the total number of hybrids over the total number of birds in counties where hybridization can occur
  hybrids_counties <- hybrid_zone_dataset %>% 
    filter(category=="hybrid")
  parental1_counties <-hybrid_zone_dataset %>% 
    filter(common_name == hybrid_system[1])
  parental2_counties <- hybrid_zone_dataset %>% 
    filter(common_name == hybrid_system[2])
  total_counties <- hybrid_zone_dataset %>% 
    filter(category=="hybrid"| category=="species")
  
  total_rate_counties <- (sum(hybrids_counties$observation_count))/(sum(total_counties$observation_count))
  parental1_rate_counties <- (sum(hybrids_counties$observation_count))/(sum(parental1_counties$observation_count))
  parental2_rate_counties <- (sum(hybrids_counties$observation_count))/(sum(parental2_counties$observation_count))
  
  
  ###3) hybridization rates in a defined geographic area
    #NOTE: this method is really only relevant for the Great Plains hybrid systems, and only works if
    #the eastern species is listed before the western species in the original character vectors
  
  #filtering down to just parental 1
  parental1 <- filtered_dataset %>%
    filter(common_name==hybrid_system[1]) %>%
    filter(observation_count >= 3)
  #finding the westernmost longitude in this filtered dataset
  west_limit <- min(parental1$longitude) - 1
  
  #filtering down to just parental 2
  parental2 <- filtered_dataset %>%
    filter(common_name==hybrid_system[2]) %>%
    filter(observation_count >= 3)
  #finding the easternmost longitude in this filtered dataset
  east_limit <- max(parental2$longitude) + 1
  
  #restricting the original dataset to be hemmed in by this longitude
  hybrid_area <- filtered_dataset %>%
    filter(longitude > west_limit & longitude < east_limit)
  
  hybrids_geographic <- hybrid_area %>% 
    filter(category=="hybrid")
  parental1_geographic <- hybrid_area %>% 
    filter(common_name==hybrid_system[1])
  parental2_geographic <- hybrid_area %>% 
    filter(common_name==hybrid_system[2])
  total_geographic <- hybrid_area %>% 
    filter(category=="hybrid"| category=="species")
  
  total_rate_geographic = (sum(hybrids_geographic$observation_count))/(sum(total_geographic$observation_count))
  parental1_rate_geographic = (sum(hybrids_geographic$observation_count))/(sum(parental1_geographic$observation_count))
  parental2_rate_geographic = (sum(hybrids_geographic$observation_count))/(sum(parental2_geographic$observation_count))
  
  
  ###bringing everything together in one dataframe
  names <- c("all", hybrid_system[1], hybrid_system[2])
  rangewide <- c(total_rate_rangewide, parental1_rate_rangewide, parental2_rate_rangewide)
  counties <- c(total_rate_counties, parental1_rate_counties, parental2_rate_counties)
  geographic <- c(total_rate_geographic, parental1_rate_geographic, parental2_rate_geographic)
  
  hybrid.data <- data.frame(names, rangewide, counties, geographic, stringsAsFactors=FALSE)
  print(hybrid.data)
  }
  
  print(substitute())
}
