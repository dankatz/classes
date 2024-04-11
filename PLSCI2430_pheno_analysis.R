#preparing data for PLSCI 2430 phenology lab analysis


### set up working environment #######################################################################
library(dplyr)
library(tidyr)
library(lubridate)
library(rnpn)
library(readr)
library(here) #install.packages("Rtools")
library(purrr)
library(ggplot2)
library(googlesheets4)
library(sf)
library(prism)

#rm(list=ls())

#prism_set_dl_dir("~/prism") 
prism_set_dl_dir("C:/Users/dsk273/Documents/prism")
#prism_set_dl_dir("C:/Users/danka/Documents/prism")

setwd("C:/Users/dsk273/Box/2430_Plant Ecology and Evolution (Chelsea Specht)/Ecology/Phenology lab") 
#setwd("C:/Users/danka/Box/2430_Plant Ecology and Evolution (Chelsea Specht)/Ecology/Phenology lab") 

### herbarium data: add in and process ################################################################
herb_raw <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1QHGdkFVKz0HLat1AIl2sYYE3dHnVZeKN1QB3X-yxbi0/edit?usp=sharing", .name_repair = "universal")
herb <- herb_raw %>% 
  mutate(date_c = lubridate::ymd(date.collected..mm.dd.yyyy.),
         open_flow = case_when(open.flowers...y.n. == "Y" ~ "Y",
                               open.flowers...y.n. == "y" ~ "Y",
                               open.flowers...y.n. == "N" ~ "N",
                               open.flowers...y.n. == "n" ~ "N",
                               FALSE ~ "uh oh"),
         DOY = yday(date_c),
         year_c = year(date_c)) %>% 
  filter(!is.na(date_c)) %>% 
  filter(!is.na(open_flow))


## add in CU weather station data temperature data from 1893 ==================================================================
ith_met_raw <- read_csv("Ithaca_historical_met_data_daily.csv") 
ith_met <- ith_met_raw %>% 
  mutate(date_c_ymd = ymd(Date),
         date_c_mdy = mdy(Date),
         date_c = case_when(!is.na(date_c_ymd) ~ date_c_ymd,
                            !is.na(date_c_mdy) ~ date_c_mdy)) %>% 
  rename(Tmean = AvgTemperature) %>% 
  dplyr::select(date_c, Tmean) %>% 
  mutate(month_c = month(date_c),
         year_c = year(date_c),
         site = case_when(  date_c < ymd("1898/1/1") ~ "Lincoln Hall",
                            date_c >= ymd("1898/1/1") & date_c < ymd("1930/10/1") ~ "Bailey Hall",
                            date_c >= ymd("1930/10/1") & date_c < ymd("1943/6/1") ~ "Roberts Hall",
                            date_c >= ymd("1943/6/1") & date_c < ymd("1969/6/1") ~ "Caldwell Field",
                            date_c >= ymd("1969/6/1") ~ "Game Farm Rd"
                           ))

ith_met_month <- ith_met %>% 
 # filter(month_c == 3 | month_c == 4) %>% 
  filter(month_c == 3 ) %>% 
  group_by(year_c) %>% 
  summarize(Tmean_mar = round(mean(Tmean, na.rm = TRUE), 2))

ggplot(ith_met_month, aes(x = year_c, y = Tmean_mar)) + geom_point() +
  geom_line(aes(x = year_c, y = zoo::rollmean(Tmean_mar, 20, na.pad = TRUE)))

ith_met %>% 
  #filter(month_c == 3 ) %>% 
  group_by(year_c, site) %>% 
  summarize(Tmean_mar = mean(Tmean, na.rm = TRUE)) %>% 
  filter(!year_c %in% c(1898, 1930, 1943, 1969, 2023)) %>% 
ggplot(aes(x = year_c, y = Tmean_mar, color = site)) + geom_point() +
  geom_line(aes(x = year_c, y = zoo::rollmean(Tmean_mar, 5, na.pad = TRUE))) + theme_bw(base_size = 26) + xlab("year") + ylab("Mean temperature (°F)") 

ith_met %>% 
  #filter(month_c == 3 ) %>% 
  group_by(year_c, site) %>% 
  summarize(Tmean_mar = mean(Tmean, na.rm = TRUE)) %>% 
  filter(year_c > 1969 & year_c < 2023) %>% 
  #filter(!year_c %in% c(1898, 1930, 1943, 1969, 2023)) %>% 
  ggplot(aes(x = year_c, y = Tmean_mar, color = site)) + geom_point() +
  geom_line(aes(x = year_c, y = zoo::rollmean(Tmean_mar, 5, na.pad = TRUE))) + theme_bw(base_size = 26) + xlab("year") + ylab("Mean temperature (°F)") +
  geom_smooth(method = "lm")

test <- ith_met %>% 
  #filter(month_c == 3 ) %>% 
  group_by(year_c, site) %>% 
  summarize(Tmean_mar = mean(Tmean, na.rm = TRUE)) %>% 
  filter(year_c > 1969 & year_c < 2023) 

summary(lm(Tmean_mar ~year_c, data = test))

## download CRUTS =====================================================================================================
# https://www.nature.com/articles/s41597-020-0453-3
# using v 4.07
# install.packages("cruts")
library(cruts);library(sp);library(raster);library(stringr);library(lubridate);library(ncdf);library(rgdal)

#location of cru ts file downloaded (from https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.07/cruts.2304141047.v4.07/tmp/)
fln <- "C:/Users/danka/Box/2430_Plant Ecology and Evolution (Chelsea Specht)/Ecology/Phenology lab/cru_ts4.07.1901.2022.tmp.dat.nc"

#convert the whole netcdf to raster; this is inefficient, but the other pieces of the package don't seem to work
# and I don't want to bother messing with netcdf4
rasta <- cruts2raster(fln, timeRange=c("1901-01-01","2023-01-01"))

ith_sp <- SpatialPoints(coords = cbind(-76.47, 42.44))
crs(ith_sp) <- crs(rasta)

cruts_tmean_raw <- raster::extract(rasta, ith_sp)
name_vector <- unlist(labels(cruts_tmean_raw))[-1]
cruts_tmean_df <- data.frame(date_full = name_vector, cruts_tmean = as.vector(cruts_tmean_raw)) %>% 
  mutate(date_c_no_x = gsub("X", "", date_full),
         date_c = ymd(date_c_no_x),
         year_c = year(date_c),
         month_c = month(date_c))

cruts_tmean_df_mar_apr <- cruts_tmean_df %>% 
  filter(month_c == 3 | month_c == 4) %>% 
  group_by(year_c) %>% 
  summarize(Tmean_crut_mar_apr = mean(cruts_tmean, na.rm = TRUE))


ggplot(cruts_tmean_df_mar_apr, aes(x = year_c, y = Tmean_crut_mar_apr)) + geom_point() + geom_smooth(method = "lm") + theme_bw(base_size = 24) + xlab("year") + ylab("mean temperature in March and April (°C)")
fit <- lm(cruts_tmean_df_mar_apr$Tmean_crut_mar_apr ~ cruts_tmean_df_mar_apr$year_c)
summary(fit)

## add in weather station and CRU TS temperature data ==========================================================================
#add temperature data to herbarium data
herb <- left_join(herb, ith_met_month)
herb <- left_join(herb, cruts_tmean_df_mar_apr) 


#data vis and regression
#compare the two temperature time series
ggplot(herb, aes(x = Tmean_crut_mar_apr * (9/5) + 32, y = Tmean_mar , color = year_c)) + geom_point(size = 4) + geom_abline(slope = 1, intercept = 0, lty = 2) + xlab("Modeled temperature (°F)") + ylab("Measured temperature (°F)") + theme_bw(base_size = 24)


ggplot(herb, aes(x = Tmean_mar_apr, y = DOY, color = open_flow)) + geom_point() + facet_wrap(~Species) + geom_smooth(method = "lm")

ggplot(herb, aes(x = Tmean_crut_mar_apr, y = DOY, color = open_flow)) + geom_point() + facet_wrap(~Species) + geom_smooth(method = "lm")


ggplot(herb, aes(x = DOY, y = open_flow)) + geom_point() + theme_bw() + facet_wrap(~Species)
ggplot(herb, aes(x = year_c, y = DOY, color = open_flow)) + geom_point() + theme_bw() + facet_wrap(~Species)

ggplot(herb, aes(x = year_c, y = how.many.open.flowers..n.)) + geom_point() + theme_bw() + facet_wrap(~Species)+ geom_smooth(method = "lm")


### NPN obs: add in and process #######################################################################
#npn_plants <- npn_species(kingdom = "Plantae")

#npn_direct <- read_csv( here("data", "NPN_220620.csv")) #try reading in data if it's already downloaded
#if(exists("npn_direct") == FALSE){ #does the npn_direct object exist? If not, download it:
  npn_direct <- npn_download_status_data(
    request_source = 'Daniel Katz, Cornell',
    species_ids = c(6, #Caltha palustris
                  161, #Erythronium americanum
                  849, #Dicentra canadensis
                  850, #Dicentra cucullaria
                  1175, #Cardamine concatenata
                  2114), #Cardamine diphylla
    years = c(as.character(2009:2023)), #years to include
    phenophase_ids = c(501), #angiosperms: 501 == "Open flowers", 502 == "Pollen release (flowers)" #conifers: 495 ==  503 ==
    additional_fields = c("Observed_Status_Conflict_Flag", "partner_group")
  )
#  write_csv(npn_direct, here("data", "NPN_220620.csv"))
#}

#how many obs did each group collect
#test <- 
  npn_direct %>% 
  mutate(ithaca = case_when(latitude > 42 & latitude < 43 &
                              longitude > -77 & longitude < -76 ~ "Ithaca",
                            TRUE ~ "other"),
         yr = year(observation_date)) %>% 
  filter(ithaca == "Ithaca" & yr == 2023)  %>% 
  group_by(genus, species) %>% 
  summarize(n = n())




# flowering intensity value (which is entered about 82% of the time)
npn_active_flow <- npn_direct %>%
  #filter(phenophase_ids == 501) %>% 
  filter(phenophase_status != -1) %>% #removing observations where the observer was unsure whether the phenophase was occurring
  filter(intensity_value != -9999) %>% 
  mutate(
    yr = year(observation_date),
    flow_prop = case_when(
      phenophase_status == 0  ~ 0,
      intensity_value == "Less than 5%" ~ 0.025,
      intensity_value ==  "5-24%"~ (0.05+0.24)/2,
      intensity_value == "25-49%" ~ (0.25+0.49)/2,
      intensity_value == "50-74%" ~ (0.5+0.74)/2,
      intensity_value == "75-94%" ~ (0.75+0.94)/2,
      intensity_value == "95% or more" ~ 0.97))
#hist(npn_active_flow$flow_prop)




### average temperature of March - April ------------------------------------------------

## download monthly prism data
prism_set_dl_dir("C:/Users/danka/Documents/prism")

#get_prism_dailys(type = "tmean", minDate = "2009-01-01",  maxDate = "2023-04-30")
get_prism_monthlys(type = "tmean", years = 2009:2023, mon = 3:4, keepZip = TRUE)

prism_archive_subset("tmean", "monthly", years = 2009:2023, mon = 3:4)
#prism_archive_subset("tmean", "daily", years = 2010, mon = 3:4)



npn_active_flow_yr_sf <- npn_active_flow %>% 
  dplyr::select(longitude, latitude, site_id, yr) %>% 
  distinct() %>% 
  st_as_sf(coords = c( "longitude", "latitude"), crs = 4326) 

#prism_set_dl_dir("~/prism")
#get_prism_monthlys(type = "tmean", years = 2009:2021, keepZip = TRUE, mon = 1:4)

#function to extract time period for each year of data
extract_temp_yr <- function(focal_yr){ #focal_yr <- 2009
  tmean_rast_yr_mo <- prism_archive_subset(temp_period = "monthly", type = "tmean", years = focal_yr, mon = 3:4)
  tmean_rast2_yr_mo <- pd_stack(tmean_rast_yr_mo)
  r_mean <- raster::calc(tmean_rast2_yr_mo, mean) #raster::plot(r_mean)
  
  npn_active_flow_yr_sf_focal_yr <- filter(npn_active_flow_yr_sf, yr == focal_yr)
  tmean_data <- unlist(raster::extract(x = r_mean, 
                                       y = npn_active_flow_yr_sf_focal_yr)) %>% as.data.frame()  
  
  npn_active_flow_yr_sf_focal_yr2 <- npn_active_flow_yr_sf_focal_yr %>%  
    dplyr::select(site_id, yr) %>% 
    mutate(tmean_mar_apr = as.numeric(unlist(tmean_data)))
  
  npn_active_flow_yr_sf_focal_yr2$geometry <- NULL
  print(focal_yr)
  return(npn_active_flow_yr_sf_focal_yr2)
}

#extract_temp_yr(focal_yr = 2010, mo_start = 1, mo_end = 4)
spring_temp_site_yr <- purrr::map_dfr(.x = 2011:2023, .f = extract_temp_yr)
npn_active_flow <- left_join(npn_active_flow, spring_temp_site_yr)

npn_active_flow %>% filter(phenophase_status == 1) %>% 
  mutate(ithaca = case_when(latitude > 42 & latitude < 43 &
                          longitude > -77 & longitude < -76 ~ "Ithaca",
                          TRUE ~ "other")) %>% 
  filter(day_of_year < 200) %>% 
ggplot(aes(x = tmean_mar_apr, y = day_of_year)) + geom_point(aes(color = flow_prop)) + theme_bw() +
  facet_wrap(~common_name) + geom_smooth(method = "lm")





### Mundy obs: add in and process ######################################################################
library(dplyr)
library(tidyr)
library(lubridate)
mundy_raw <- read_csv("mundy_obs_long3.csv")

mundy <- mundy_raw %>% 
    mutate(variable = rep(c("date_start", "date_stop"), nrow(mundy_raw) / 2), 
           key = rep(1:(nrow(mundy_raw) / 2), each = 2)) %>%
    pivot_wider(id_cols = c(key,species, common_name), names_from = variable, values_from = date_start_stop) %>% 
    dplyr::select(-key) %>% 
  mutate(date_start = mdy(date_start),
         date_stop = mdy(date_stop),
         date_mean = (date_stop - date_start)/2 + date_start,
         year_c = year(date_mean),
         doy_mean = yday(date_mean),
         doy_start = yday(date_start),
         doy_stop = yday(date_stop)) 

mundy <- left_join(mundy, ith_met_month) %>%
  filter(species != "Cardamine maxima")

ggplot(mundy, aes(x = year_c, y = doy_stop)) + geom_point() + theme_bw() + geom_smooth(method = "lm", se = FALSE) +facet_wrap(~species)

ggplot(mundy, aes(x = Tmean_mar, y = doy_start)) + geom_point() + theme_bw() + geom_smooth(method = "lm") +facet_wrap(~species) +
  geom_vline(xintercept = 37.2, color = "red") + xlab("average temperate in March (F)") + ylab("start of flowering (Julian day)") +
  geom_hline(yintercept = 111, color = "purple")



### save all files for the students ####################################################################

##herbarium data
herb_export <- herb %>% 
  rename(date_collected = date_c,
         flowers_or_flower_buds = Flowers.or.flower.buds..y.n.,
         open_flowers = open_flow,
         how_many_open = how.many..n.,
         Tmean_CU_F_mar_apr = Tmean_mar_apr,
         Tmean_CRUT_C_mar_apr = Tmean_crut_mar_apr
         ) %>% 
  relocate(Species, date_collected, DOY, year_c, Tmean_CU_F_mar_apr, Tmean_CRUT_C_mar_apr, open_flowers) %>% 
  arrange(Species, open_flowers, date_collected)

herb_export

## npn data
npn_export <- npn_active_flow %>% filter(phenophase_status == 1) %>% 
  mutate(collected_in_ithaca = case_when(latitude > 42 & latitude < 43 &
                              longitude > -77 & longitude < -76 ~ "Ithaca",
                            TRUE ~ "other")) %>% 
  filter(day_of_year < 200) %>% 
  mutate(Species = paste(genus, species) ) %>% 
  relocate(Species, observation_date, DOY = day_of_year, year_c = yr, tmean_mar_apr, intensity_value, flow_prop) %>% 
  dplyr::select(-genus, - species)

npn_export


## mundy data
mundy_export <- mundy %>% rename(Species = species)


## export a file for each species for each dataset
setwd("C:/Users/dsk273/Box/2430_Plant Ecology and Evolution (Chelsea Specht)/Ecology/Phenology lab/for class to analyze 2024 b") 
sp_list <- unique(mundy_export$Species)

for(i in 1:6){
  # herb_export_sp <- filter(herb_export, Species == sp_list[i])
  # write_csv(herb_export_sp, 
  #           paste(sp_list[i], "herbarium data.csv"))
  # 
  # npn_export_sp <- filter(npn_export, Species == sp_list[i])
  # write_csv(npn_export_sp, 
  #           paste(sp_list[i], "NPN data.csv"))
  
  mundy_export_sp <- filter(mundy_export, Species == sp_list[i])
  write_csv(mundy_export_sp, paste(sp_list[i], "Mundy data_2024_c.csv"))
}

write_csv(mundy_export, "C:/Users/dsk273/Box/2430_Plant Ecology and Evolution (Chelsea Specht)/Ecology/Phenology lab/for class to analyze 2024 b/all_species.csv")

