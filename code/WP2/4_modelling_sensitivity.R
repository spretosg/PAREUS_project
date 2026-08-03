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
library(SSDM)


source("WP2/wp2_functions_utils.R")
stud_id<-"FRL04"


main_dir<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",stud_id,"/raw_data_backup")
eval_round<-"R1" #R2



## load and stack the predictor variables
in_var <- list.files(paste0(main_dir,"/env_var"), full.names = TRUE)
pred<-SSDM::load_var(path=paste0(main_dir,"/2_env_var"), categorical = "lulc")



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
mapper<-read.csv(paste0(main_dir,"/mapper.csv"))%>%filter(siteID == stud_id & userID %in% users)

es<-unique(ind_pols$esID)

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

results_list <- list()
varimp_list<-list()
es_ids<-es[1]
t0<-Sys.time()
for (es in es_ids) {
  
  message("Processing ES: ", es)
  
  pol_es <- ind_pols %>% filter(esID == es)

  # ---- FULL PTS ----
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

  algos <- c("GLM", "RF", "GAM", "Maxent")
  
  mods <- lapply(algos, function(a)
    SSDM::modelling(
      a,
      df_full,
      pred_w,
      Xcol = "lon",
      Ycol = "lat",
      cv = "holdout",
      cv.param = c(0.7, 2),
      final.fit.data = "all"
    )
  )
  
  names(mods) <- algos
  perf <- do.call(rbind,
                  lapply(mods, function(x) x@evaluation)
  )

  imp_norm <- lapply(mods, function(m){
    x <- m@variable.importance
    x / sum(x)
  })
  
  alg_stack <- list(
    terra::rast(mods$GLM@projection),
    terra::rast(mods$RF@projection),
    terra::rast(mods$GAM@projection)
    #terra::rast(mods$Maxent@projection)
  )
  names(alg_stack) <- c("GLM","RF","GAM")
  
  # algo_stack <- terra::rast(c(r_RF, r_GLM, r_GAM))
  # 

  terra::writeRaster(
    algo_stack,
    paste0(stud_id,"_", es, "_algoStack.tif"),
    overwrite = TRUE
  )
  
  
  results_list[[as.character(es)]] <- user_scores
  projection_results[[as.character(es)]] <- projection_list
  
  preds <- rast(projection_results[[as.character(es)]])
  
  names(preds) <- names(projection_results[[as.character(es)]])
  # 
  # mean_pred <- terra::app(preds, mean, na.rm = TRUE)
  # sd_pred <- terra::app(preds, sd, na.rm = TRUE)
  # cv_pred<-sd_pred/mean_pred
  # 
  # range_pred <- terra::app(preds, function(x)
  #   max(x, na.rm=TRUE) - min(x, na.rm=TRUE))


### user influence analysis


  users_es <- unique(pol_es$userID)
  
  ## Full model ------------------------------------------------------------
  
  m_full <- SSDM::modelling(
    "RF",
    df_full,
    pred_w,
    Xcol = "lon",
    Ycol = "lat",
    cv = "holdout",
    cv.param = c(0.7,2),
    final.fit.data = "all"
  )
  
  auc_full <- m_full@evaluation$AUC
  varimp_full <- m_full@variable.importance
  proj_full <- terra::rast(m_full@projection)
  
  
  
  
  
  
  ## Table to store results ------------------------------------------------
  
  user_scores <- data.frame(
    userID = users_es,
    AUC_full = auc_full,
    AUC_minusUser = NA,
    delta_AUC = NA,
    MeanDiff = NA,
    RMSE = NA,
    Corr = NA,
    Delta_dem = NA,
    Delta_lulc = NA,
    Delta_int = NA,
    Delta_acc = NA
  )
  
  ## Store prediction maps if desired
  projection_list <- list()
  projection_results<-list()
  
  ## Leave-one-user-out ----------------------------------------------------
  
  for (u in users_es) {
    
    message("Removing user: ", u)
    
    pol_minus <- pol_es %>%
      filter(userID != u)
    
    pts_minus <- make_presence_points(pol_minus, A_roi, all_back_pts)
    pts_minus <- st_transform(pts_minus, st_crs(pred))
    pts_minus_sp <- as(pts_minus, "Spatial")
    
    extracted2 <- raster::extract(pred_w, pts_minus_sp)
    
    df_minus <- cbind(pts_minus, extracted2)
    df_minus <- cbind(df_minus, st_coordinates(df_minus))
    colnames(df_minus)[colnames(df_minus) %in% c("X","Y")] <- c("lon","lat")
    
    df_minus <- st_drop_geometry(df_minus)
    df_minus <- na.omit(df_minus)
    df_minus <- df_minus[,c("lon","lat")]
    
    df_minus$SPECIES <- "pres"
    
    ## Fit model
    m_minus <- SSDM::modelling(
      "RF",
      df_minus,
      pred_w,
      Xcol="lon",
      Ycol="lat",
      cv="holdout",
      cv.param=c(0.7,2),
      final.fit.data="all"
    )
    
    ## Evaluation
    auc_minus <- m_minus@evaluation$AUC
    
    ## Variable importance
    varimp_minus <- m_minus@variable.importance
    
    ## Prediction map
    proj_minus <- terra::rast(m_minus@projection)
    
    projection_list[[as.character(u)]] <- proj_minus
    
    ## Spatial difference
    diff <- proj_full - proj_minus
    
    mean_diff <- terra::global(abs(diff), "mean", na.rm=TRUE)[1,1]
    
    rmse <- sqrt(
      terra::global(diff^2, "mean", na.rm=TRUE)[1,1]
    )
    
    vals <- na.omit(cbind(
      terra::values(proj_full),
      terra::values(proj_minus)
    ))
    
    corr <- cor(vals[,1], vals[,2])
    
    ## Save
    
    user_scores[user_scores$userID==u,"AUC_minusUser"] <- auc_minus
    user_scores[user_scores$userID==u,"delta_AUC"] <- auc_full-auc_minus
    
    user_scores[user_scores$userID==u,"MeanDiff"] <- mean_diff
    user_scores[user_scores$userID==u,"RMSE"] <- rmse
    user_scores[user_scores$userID==u,"Corr"] <- corr
    
    user_scores[user_scores$userID==u,"Delta_dem"] <-
      varimp_full[1]-varimp_minus[1]
    
    user_scores[user_scores$userID==u,"Delta_lulc"] <-
      varimp_full[2]-varimp_minus[2]
    
    user_scores[user_scores$userID==u,"Delta_int"] <-
      varimp_full[3]-varimp_minus[3]
    
    user_scores[user_scores$userID==u,"Delta_acc"] <-
      varimp_full[4]-varimp_minus[4]
    
  }
  
  results_list[[as.character(es)]] <- user_scores
  projection_results[[as.character(es)]] <- projection_list
  user_stack <- terra::rast(projection_list)
  
  names(user_stack) <- names(projection_list)
  
  terra::writeRaster(
    user_stack,
    filename = file.path(
      out_dir,
      paste0(stud_id,"_", es, "_userStack.tif")
    ),
    overwrite = TRUE
  )
  

  
  preds <- rast(projection_results[[as.character(es)]])
  
  names(preds) <- names(projection_results[[as.character(es)]])
  
  # mean_pred <- app(preds, mean, na.rm=TRUE)
  # sd_pred <- app(preds, sd, na.rm=TRUE)
  # 
  # 
  # range_pred <- app(preds,
  #                   function(x)
  #                     max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
  # 
  # 
  # uncertainty <- app(preds,
  #                    function(x)
  #                      sd(x)/mean(x))
  # uncertainty <- app(preds,
  #                    function(x)
  #                      quantile(x,0.75,na.rm=T)-quantile(x, 0.25, na.rm=T))

} 

