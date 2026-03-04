## post hoc analysis

## this script analyses the importance of single participant
## the importance of a single predictor
## uses declustering techniques

library(sf)
library(leaflet)
library(DT)
library(dplyr)
library(ggplot2)
library(dplyr)
library(terra)

source("WP2/wp2_functions_utils.R")
stud_id<-"FRL04"
main_dir<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",stud_id,"/raw_data_backup")
eval_round<-"R1" #R2



## load and stack the predictor variables
in_var <- list.files(paste0(main_dir,"/env_var"), full.names = TRUE)
pred<-SSDM::load_var(path=paste0(main_dir,"/env_var"))

# 2. Read all rasters
rasters <- lapply(in_var, rast)


## read participants data
ind_pols<-read_sf(paste0(main_dir,"/ind_polys_",eval_round,".gpkg"))%>%dplyr::filter(siteID == stud_id)
stud_site<-read_sf(paste0(main_dir,"/study_site.gpkg"))%>%dplyr::filter(siteID == stud_id)
A_roi<-stud_site$siteAREAkm2*10^6
resolution = 250^2
all_back_pts<- round(A_roi/resolution,0)

## join imp_access to ind_pols via esID and userID from esmappingR1.csv
es_mapping<-read.csv(paste0(main_dir,"/es_mappingR1.csv"))

ind_pols <- ind_pols %>%
  dplyr::left_join(es_mapping %>% dplyr::select(esID, userID, imp_acc),
            by = c("esID", "userID"))

users<-unique(ind_pols$userID)
es<-unique(ind_pols$esID)
#1. at es level
#1.1 run a model with all users -- var imp
#1.2. run a model with hold out for user 
#
#Sensitivity: modelling parameters, 
# min_in_pts<-10
# for(n in 1:length(es)){
#   tmp_es<-es[n]
#   polygon<-ind_pols%>%filter(esID == tmp_es)
#   for (i in 1:nrow(polygon)) {
#     A_tmp <- as.numeric(st_area(polygon[i,]))
#     tmp_ratio<-A_tmp/A_roi
#     tmp_pts<-round(all_back_pts*tmp_ratio,0)
#     
#     if(tmp_pts<=min_in_pts){
#       tmp_pts<-min_in_pts
#     }else{
#       tmp_pts<-tmp_pts
#     }
#     # npts in this poly must be max_pts*tmp_ratio*es_value
#     #+1 its 0-5 scale
#     tmp_es_val<-((1+polygon[i,]$es_value)/5)
#     tmp_pts = st_sample(polygon[i,], round(tmp_pts*tmp_es_val,0),type="random")
#     tmp_pts<-st_as_sf(tmp_pts)
#     tmp_pts$inside<-rep(1,nrow(tmp_pts))
#     if(i==1){
#       pts_in<-tmp_pts
#     }else{
#       pts_in<-rbind(pts_in,tmp_pts)
#     }
#     
#   }
#   
#   weight_access <- as.numeric(as.numeric(mean(polygon$imp_acc))/5)
#   pred_w<-raster::stack(pred$dem*1, pred$lulc*1, pred$int*1, pred$acc*weight_access)
#   pts_in<-st_transform(pts_in,st_crs(pred))
#   pts_in_sp <- as(pts_in, "Spatial")
#   
#   # 3. Extract raster values
#   extracted_values <- raster::extract(pred_w, pts_in_sp)
#   
#   # 4. Combine extracted values with original points
#   # If pred_w has multiple layers, extract returns a matrix
#   if (is.matrix(extracted_values)) {
#     colnames(extracted_values) <- names(pred_w)
#     pts_in <- cbind(pts_in, extracted_values)
#   } else {
#     # Single layer case
#     pts_in <- cbind(pts_in, raster_value = extracted_values)
#   }
#   
#   # 5. Add coordinates
#   pts_in <- cbind(pts_in, st_coordinates(pts_in))
#   colnames(pts_in)[colnames(pts_in) %in% c("X", "Y")] <- c("lon", "lat")
#   
#   # 6. Convert to data frame and remove NAs
#   pts <- st_drop_geometry(pts_in)
#   pts <- na.omit(pts)
#   pts <- pts %>% dplyr::select(lon, lat)
#   
#   #######################
#   
#   
#   # pts <- do.call(rbind, st_geometry(pts_in)) %>% 
#   #   as_tibble() %>% setNames(c("lon","lat"))
#   pts$SPECIES<-rep("pres",nrow(pts))
#   
#   
#   # Train model
#   SDM <- SSDM::modelling(
#     'RF', pts, 
#     pred_w,
#     Xcol = 'lon', Ycol = 'lat',
#     cv = "holdout",
#     cv.param = c(0.7, 2),
#     final.fit.data = "all"
#   )
#   
#   
# }


### ---- helper fct
make_presence_points <- function(polys, A_roi, all_back_pts, min_in_pts = 10) {
  
  pts_list <- lapply(1:nrow(polys), function(i) {
    
    A_tmp <- as.numeric(st_area(polys[i,]))
    prop <- A_tmp / A_roi
    npts <- round(all_back_pts * prop)
    npts <- max(npts, min_in_pts)
    
    # ES value scaling
    es_scale <- (1 + polys$es_value[i]) / 5
    npts <- round(npts * es_scale)
    
    # sample points
    pts_sf <- st_sample(polys[i,], npts, type = "random") |> 
      st_as_sf() |> 
      mutate(inside = 1)
    
    pts_sf
  })
  
  do.call(rbind, pts_list)
}


### ---- leave one out analysis
library(SSDM)

results_list <- list()
varimp_list<-list()
es_ids<-es[1]
t0<-Sys.time()
for (es in es_ids) {
  
  message("Processing ES: ", es)
  
  pol_es <- ind_pols %>% filter(esID == es)

  
  # ---- FULL MODEL ----
  pts_full <- make_presence_points(pol_es, A_roi, all_back_pts, min_in_pts = 10)
  pts_full <- st_transform(pts_full, st_crs(pred))
  pts_full_sp <- as(pts_full, "Spatial")
  
  pred_w <- raster::stack(
    pred$dem*1,
    pred$lulc*1,
    pred$int*1,
    pred$acc * (mean(pol_es$imp_acc)/5)
  )
  
  extracted <- extract(pred_w, pts_full_sp)
  df_full <- cbind(pts_full, extracted)
  
  df_full <- cbind(df_full, st_coordinates(df_full))
  colnames(df_full)[colnames(df_full) %in% c("X", "Y")] <- c("lon", "lat")
  
  # 6. Convert to data frame and remove NAs
  df_full <- st_drop_geometry(df_full)
  df_full <- na.omit(df_full)
  df_full <- df_full %>% dplyr::select(lon, lat)
  

  df_full$SPECIES <- "pres"
  
  m_full <- SSDM::modelling("RF", df_full,
                            pred_w, Xcol="lon", Ycol="lat",
                            cv="holdout", cv.param=c(0.7,2),
                            final.fit.data="all")
  
  auc_full <- m_full@evaluation$AUC
  varimp_es<-m_full@variable.importance
  varimp_list[[as.character(es)]] <- varimp_es
  
  rf_pred_rast <- terra::rast(m_full@projection)
  #write_modelled raster to drive
  raster_out<-terra::writeRaster(rf_pred_rast,paste0(main_dir,"/4_mean_R1/",es,"_RF_all.tif"))
  
  # ## compute corr between RF_all and mean
  # mean_rast <- terra::rast(paste0(main_dir, "/4_mean_R1/", es, "_mean.tif"))
  # 
  # # Ensure same resolution, extent, CRS
  # mean_rast <- terra::resample(mean_rast, rf_pred_rast)
  # 
  # # Extract values
  # val_rf   <- terra::values(rf_pred_rast, mat = FALSE)
  # val_mean <- terra::values(mean_rast, mat = FALSE)
  # 
  # # Remove missing pairs
  # keep <- complete.cases(val_rf, val_mean)
  # 
  # corr_full <- cor(val_mean[keep], val_rf[keep])
  # ---- USER REMOVAL ----
  users_es <- unique(pol_es$userID)
  
  user_scores <- data.frame(
    userID = users_es,
    AUC_full = auc_full,
    AUC_minusUser = NA,
    delta_AUC = NA
  )
  
  for (u in users_es) {
    message("   Removing user: ", u)
    
    pol_minus <- pol_es %>% filter(userID != u)

    pts_minus <- make_presence_points(pol_minus, A_roi, all_back_pts)
    pts_minus <- st_transform(pts_minus, st_crs(pred))
    pts_minus_sp <- as(pts_minus, "Spatial")
    
    extracted2 <- extract(pred_w, pts_minus_sp)
    df_minus <- cbind(pts_minus, extracted2)
    df_minus <- cbind(df_minus, st_coordinates(df_minus))
    colnames(df_minus)[colnames(df_minus) %in% c("X", "Y")] <- c("lon", "lat")
    
    # 6. Convert to data frame and remove NAs
    df_minus <- st_drop_geometry(df_minus)
    df_minus <- na.omit(df_minus)
    df_minus <- df_minus %>% dplyr::select(lon, lat)

    df_minus$SPECIES <- "pres"
    
    m_minus <- SSDM::modelling("RF", df_minus,
                               pred_w, Xcol="lon", Ycol="lat",
                               cv="holdout", cv.param=c(0.7,2),
                               final.fit.data="all")
    varimp_user<-m_minus@variable.importance
    
    auc_minus <- m_minus@evaluation$AUC
    
    user_scores[user_scores$userID == u, "AUC_minusUser"] <- auc_minus
    user_scores[user_scores$userID == u, "delta_AUC"]      <- auc_full - auc_minus
    user_scores[user_scores$userID == u, "varimp_dem"]  <- varimp_user[1]
    user_scores[user_scores$userID == u, "varimp_lulc"]  <- varimp_user[2]
    user_scores[user_scores$userID == u, "varimp_int"]  <- varimp_user[3]
    user_scores[user_scores$userID == u, "varimp_acc"]  <- varimp_user[4]
  }
  
  results_list[[as.character(es)]] <- user_scores
}
print(Sys.time()-t0)
# weight predictors

# pred_w<-stack(pred$dem*1, pred$eii*1, pred$acc*as.numeric(input$imp_acc))




