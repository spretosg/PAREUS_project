## test data read tif
## make pu of stud area
## sample cost, es_cond and es_capacity per PU

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
source("code/WP4/wp4_functions_utils.R")
target_site<-"SK021"


####---- Input and processing ----####
## read PU from 01_pa_status.R
grid<-st_read(paste0("outputs/WP4/01_PA_analysis/",target_site,"_input_grid.json"))

### Features
## Mean ecosystem service (cookbook nr.2)
grand_mean_es<-terra::rast(paste0("data/WP2/",target_site,"_Wmean.tif"))
# ES groups
es_raster <- rast(paste0("data/WP4/es_individual_",target_site,".tif"))

es_groups <- list(
  provisioning = c("wild_plant", "wild_hunt", "farm", "mat"),
  regulating   = c("erosion", "flood", "habitat"),
  cultural     = c("aest", "recr","sense")
)

#prov es mean
mean_prov <- mean(es_raster[[es_groups$provisioning]])
mean_reg  <- mean(es_raster[[es_groups$regulating]])
mean_cult <- mean(es_raster[[es_groups$regulating]])

## Ecosystem condition (cookbook nr.2)
es_cond<-terra::rast(paste0("data/WP2/",target_site,"_ECSI.tif"))

### costs
##based on ecosystem service exclusiveness
cost_es<-terra::rast(paste0("data/WP4/",target_site,"_cost_raster_es.tif"))
## based on policy coherence (cookbook nr. 3)
cost_policy_file <- file.path(
  "data/WP4/",
  paste0(target_site, "_cost_raster_pol.tif")
)

cost_policy <- if (file.exists(cost_policy_file)) {
  terra::rast(cost_policy_file)
} else {
  NULL
}


target_crs<-st_crs(grid)$wkt

rasters <- list(
  cost_es       = cost_es,
  cost_policy   = cost_policy,
  grand_mean_es = grand_mean_es,
  es_cult       = mean_cult,
  es_reg        = mean_reg,
  es_prov       = mean_prov,
  eco_cond = es_cond
)

rasters <- Filter(Negate(is.null), rasters)

list2env(rasters, envir = environment())

grid[paste0("sampled_", names(rasters))] <-
  lapply(rasters, \(r)
         terra::extract(r, grid, fun = mean, na.rm = TRUE)[, 2]
  )

#calculate distance from cells outside core PA to core PA as defined in step 01_pa_status.R
out   <- grid%>%filter(existing_corePA == F)
pa_core <- grid%>%filter(existing_corePA == T)

out$min_distance <- apply(st_distance(out, pa_core), 1, min)
pa_core$min_distance<-0
grid<-rbind(out,pa_core)

## remove na grid cells based on features and costs (not)
cols_to_check <- c(
  "sampled_grand_mean_es",
  "sampled_es_cult",
  "sampled_es_prov",
  "sampled_es_reg",
  "sampled_cost_es",
  "sampled_eco_cond"
)

if ("sampled_cost_policy" %in% names(grid)) {
  cols_to_check <- c(cols_to_check, "sampled_cost_policy")
}

grid_clean <- grid %>%
  filter(!if_any(all_of(cols_to_check), ~ is.na(.) | is.nan(.)))

## scale
if ("sampled_cost_policy" %in% names(grid_clean)) {
  cols_to_check <- c(cols_to_check, "sampled_cost_policy")
}

grid_clean <- zero_one_scale(
  grid_clean,
  cols = cols_to_check
)


st_write(grid_clean, paste0("outputs/WP4/02_optim/",target_site,"_input_final_grid.json"), driver = "GeoJSON", overwrite = T)
