### calculate mean raster per ES

library(terra)
library(sf)
library(dplyr)
library(spatstat.geom)
library(spatstat.model)
source("WP2/wp2_functions_utils.R")
stud_id<-"FRL04"
main_dir<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",stud_id,"/raw_data_backup")
eval_round<-"R1" #R2

if(eval_round == "R1"){
  target_dir<-"/3_ind_R1"
}else{
  target_dir <-"/5_ind_R2"
}


es_ratings<-read.csv(paste0(main_dir,"/ahp_weights.csv"))

# List the first 10 subfolders
subfolders <- list.dirs(paste0(main_dir,target_dir), full.names = TRUE, recursive = FALSE)[1:10]

## ind polygons to create kernel density surface
ind_pols<-read_sf(paste0(main_dir,"/ind_polys_",eval_round,".gpkg"))%>%dplyr::filter(siteID == stud_id)
#polys_proj <- st_transform(ind_pols, 3857)



#### ---- create a kernel density map for each es ####
# 1. Create centroids and a numeric weight
# polys_proj<-ind_pols%>%filter(esID =="recr")
# pts <- st_centroid(polys_proj)
# 
# 
# ### IDW
# # Create interpolation grid (choose resolution)
# r <- rast(ext(mean_raster), res = 0.005) 
# crs(r) <- st_crs(pts)$wkt
# xy <- terra::xyFromCell(r, 1:ncell(r))
# coop <- st_as_sf(as.data.frame(xy), coords = c("x", "y"),
#                  crs = st_crs(pts))
# 
# 
# ### distance based only
# # If pts has 0 rows → return pts unchanged with dist column
# if (nrow(coop) == 0) {
#   coop$mindist <- numeric(0)
# } else {
#   distmat <- st_distance(coop, pts)
#   coop$mindist <- apply(distmat, 1, min)
# }
# pred <- terra::rasterize(coop, r, field = "mindist", fun = "mean")
# # mn <- global(pred, "min", na.rm = TRUE)[1]
# # mx <- global(pred, "max", na.rm = TRUE)[1]
# pred_log <- log(pred)
# 
# # pred <- (pred - as.numeric(mn)+1) / (as.numeric(mx)-as.numeric(mn))
# pred_sum <- global(pred_log, "sum", na.rm = TRUE)[1]
# w <- pred_log / as.numeric(global(pred_log, "max", na.rm = TRUE)[1])
# w<-terra::resample(w,mean_raster,"bilinear")

rast_weights<-read.csv(paste0(main_dir,"/es_mapping",eval_round,".csv"))%>%filter(siteID == stud_id)
#### ---- Loop through each subfolder#### (we need to weight it with uncertainty!)

grand_mean<-rast()
for (folder in subfolders) {
  
  # List all raster files (you can adjust the pattern if needed)
  raster_files <- list.files(folder, full.names = TRUE)
  es_tmp<-sub(".*/", "", folder)
  rasters <- lapply(raster_files, rast)
  rasters <- terra::rast(rasters)

  ind_weights<-rast_weights%>%filter(esID==es_tmp)%>%select(confidence,userID)%>%mutate(conf_adj = case_when(confidence>0 ~confidence/5,
                                                                                                             confidence == 0 ~ 0))
  
  ind_weights <- ind_weights %>%
    filter(userID %in% names(rasters)) %>%
    mutate(idx = match(userID, names(rasters))) %>%
    arrange(idx)

  ### validate only use not 0 rasters
  valid_rasters <- lapply(raster_files, function(f) {
    r <- rast(f)
    
    # compute the maximum pixel value
    mx <- global(r, fun = "max", na.rm = TRUE)[1,1]
    
    # if max == 0 → raster contains only zeros → skip
    if (mx == 0 || is.na(mx)) {
      message("Skipping empty raster: ", f)
      return(NULL)
    }
    
    return(r)
  })
  
  # drop NULLs
  valid_rasters <- Filter(Negate(is.null), valid_rasters)
  
  # Optionally combine into a SpatRaster stack
  if (length(valid_rasters) > 0) {
    raster_stack <- rast(valid_rasters)
  }
  
  raster_stack <- raster_stack[[ind_weights$idx]]   # reorder raster too (extra safety)
  
  #confidence weighted raster (not es importance weighted!)
  r_weighted <- raster_stack * ind_weights$conf_adj
  

  # Calculate mean raster

  mean_raster <- mean(r_weighted, na.rm = TRUE)
  #weight with ahp ratings
  tmp_rating<-as.numeric(es_ratings%>%filter(esID==es_tmp)%>%select(pref_adj))
  mean_rast_w<-mean_raster*tmp_rating
  grand_mean<-c(grand_mean,mean_rast_w)
  
  crs(mean_raster) <- "EPSG:4326"
  
  cv_raster <- cv_rast(raster_stack)
  
  ## calc uncertainty weighted ES mean
  ES_benefit<-1/cv_raster*mean_raster

  # Create output filename
  # folder_name <- basename(folder)
  mean_path<-paste0(main_dir,"/4_mean_",eval_round)
  cv_path<-paste0(main_dir,"/5_cv_",eval_round)
  ben_path<-paste0(main_dir,"/6_benefit_",eval_round)
  if(!exists(mean_path)){
    dir.create(mean_path)
  }
  if(!exists(cv_path)){
    dir.create(cv_path)
  }
  if(!exists(ben_path)){
    dir.create(ben_path)
  }
  out_mean<- file.path(mean_path, paste0(basename(folder), "_mean.tif"))
  out_cv<- file.path(cv_path, paste0(basename(folder), "_cv.tif"))
  out_ben<- file.path(ben_path, paste0(basename(folder), "_ben.tif"))
  names(mean_raster)<-basename(folder)
  names(cv_raster)<-basename(folder)
  names(ES_benefit)<-basename(folder)

  writeRaster(mean_raster, out_mean, overwrite = TRUE)
  writeRaster(cv_raster, out_cv, overwrite = TRUE)
  writeRaster(ES_benefit, out_ben, overwrite = TRUE)

  cat("Saved to:", main_dir, "\n")
}

w_mean<-sum(grand_mean,na.rm=T)

names(w_mean)<-"es_weighted_mean"

writeRaster(w_mean, "P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/FRL04/es_weight_mean.tif", overwrite = TRUE)

## subset ind_polys
write_sf(ind_pols,paste0(main_dir,"/ind_polys_R1.gpkg"))
