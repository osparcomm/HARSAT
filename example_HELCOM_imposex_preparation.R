# get imposex data 

library(dplyr)
library(tibble)
library(tidyr)

wk_timeSeries <- biota_timeSeries

wk_data <- wk_timeSeries$data %>% 
  filter(
    ctsm_get_info("determinand", determinand, "group", "biota", sep = "_") %in% 
      "Imposex"
  ) %>% 
  select(c("seriesID", "station_code", "year", "species", "determinand", 
           "concentration", "noinp", "%FEMALEPOP"))

wk_data <- left_join(
  wk_data, 
  wk_timeSeries$stations[c("station_code", "country", "HELCOM_subbasin")],
  by = "station_code"
)

rm(wk_timeSeries)


# see where we have individual data 

stopifnot(round(wk_data$noinp) == wk_data$noinp)

wk_data <- wk_data %>%
  mutate(noinp = as.integer(noinp)) %>%
  group_by(station_code, determinand, species, year) %>% 
  filter(all(noinp == 1L)) %>% 
  ungroup() %>% 
  droplevels() %>% 
  select(-noinp)
  


# get cut points by country 
# have used terminology from MIME code (
# e.g. regionID = country here, whereas regionID ~ country + region for 
# MIME)

wk_data <- unite(wk_data, "regionID", country, sep = " ", remove = FALSE)


# have difficulty dealing with Neptunea because nearly all zero and gap in values 
# (no records for class 2 or 3)
# for example: 
# Denmark 28 1 0 0 1 0 0 

with(wk_data, table(regionID, concentration, species))


wk_data <- wk_data %>% 
  filter(!.data$species %in% "Neptunea antiqua") %>% 
  select(-determinand, -"%FEMALEPOP") %>%
  droplevels()

summary(wk_data)


# construct new variable VDS so we can combine larger categories with few observations without 
# changing the raw data

wk_data <- mutate(wk_data, VDS = .data$concentration)

with(wk_data, table(regionID, VDS, species))


# combine larger categories with few observations

wk_data <- within(wk_data, {
  
  id <- species == "Buccinum undatum"
  VDS[id & VDS >= 2] <- 2
  
})

with(wk_data, table(regionID, VDS, species))


# redefine indexID and regionID so don't have to worry about multiple species

wk_data <- wk_data %>% 
  unite("indexID", station_code, year, species, sep = " ", remove = FALSE) %>% 
  unite("regionID", regionID, species, sep = " ", remove = FALSE) %>% 
  mutate(indexID = factor(indexID)) %>% 
  as.data.frame()

with(wk_data, table(regionID, VDS))


# split by regionID, estimate cut points and level for each indexID

wk_split <- split(wk_data, wk_data$regionID)


library("parallel")
library("pbapply")

wk.cores <- detectCores()
wk.cluster <- makeCluster(wk.cores - 1)

clusterExport(wk.cluster, c("wk_split", ctsm.VDS.varlist))

clusterEvalQ(wk.cluster, {
  library("MASS")
})  

biota.VDS.estimates <- pblapply(
  wk_split, 
  ctsm.VDS.index.opt, 
  calc.vcov = TRUE, 
  cl = wk.cluster
)

# check convergence

all(sapply(biota.VDS.estimates, "[[", "convergence") == 0)

# saveRDS(
#   biota.VDS.estimates,
#   file.path("RData", "VDS estimates.rds")
# )


# get confidence limits on estimated VDSI for each indexID

clusterExport(wk.cluster, "biota.VDS.estimates")

biota.VDS.cl <- parLapply(wk.cluster, biota.VDS.estimates, ctsm.VDS.cl)
biota.VDS.cl <- do.call(rbind, biota.VDS.cl)

row.names(biota.VDS.cl) <- do.call(
  paste, 
  biota.VDS.cl[c("station_code", "year", "species")]
)
summary(biota.VDS.cl)

# saveRDS(
#   biota.VDS.cl,
#   file.path("RData", "VDS confidence limits.rds")
# )

stopCluster(wk.cluster)

rm(wk_data, wk_split, wk.cluster, wk.cores)
