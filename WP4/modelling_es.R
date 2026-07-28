

main_dir<-"P:/312204_pareus/"

## function to sample presence points
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


run_algorithm_comparison <- function(es_ids,
                                     studyID,
                                     ind_pols,
                                     pred,
                                     A_roi,
                                     all_back_pts,
                                     out_dir,
                                     algos = "RF") {
  
  
  performance_all <- list()
  imp_all<-list()
  
  for (es in es_ids) {
    
    message("Processing ES ", es)
    
    pol_es <- dplyr::filter(ind_pols, esID == es)
    
    if (nrow(pol_es) == 0)
      next
    pol_es<-st_make_valid(pol_es)
    ## Presence points
    pts_full <- make_presence_points(
      pol_es,
      A_roi,
      all_back_pts,
      min_in_pts = 10
    )
    
    pts_full <- sf::st_transform(pts_full, sf::st_crs(pred))
    pts_full_sp <- as(pts_full, "Spatial")
    
    ## Predictor stack
    pred_w <- raster::stack(
      pred$dem * 1,
      pred$lulc * 1,
      pred$int * 1,
      pred$acc * (mean(pol_es$imp_acc, na.rm = TRUE) / 5)
    )
    
    extracted <- raster::extract(pred_w, pts_full_sp)
    
    df_full <- cbind(pts_full, extracted)
    df_full <- cbind(df_full, sf::st_coordinates(df_full))
    
    colnames(df_full)[colnames(df_full) %in% c("X","Y")] <- c("lon","lat")
    
    df_full <- sf::st_drop_geometry(df_full)
    df_full <- na.omit(df_full)
    df_full <- dplyr::select(df_full, lon, lat)
    
    df_full$SPECIES <- "pres"
    
    ## Fit models
    mods <- lapply(algos, function(a)
      SSDM::modelling(
        a,
        df_full,
        pred_w,
        Xcol = "lon",
        Ycol = "lat",
        cv = "holdout",
        cv.param = c(0.7,2),
        final.fit.data = "all"
      )
    )
    
    names(mods) <- algos
    
    ## Performance table
    perf <- do.call(rbind,
                    lapply(names(mods), function(a){
                      
                      p <- as.data.frame(mods[[a]]@evaluation)
                      
                      p$Algorithm <- a
                      p$studyID <- studyID
                      p$ES <- es
                      
                      p
                      
                    })
    )
    
    performance_all[[as.character(es)]] <- perf
    
    ## Prediction stack
    alg_stack <- terra::rast(
      lapply(mods, function(m)
        terra::rast(m@projection))
    )
    
    names(alg_stack) <- names(mods)
    
    terra::writeRaster(
      alg_stack,
      filename = file.path(
        out_dir,
        paste0(studyID, "_", es, ".tif")
      ),
      overwrite = TRUE
    )
    
    ## Variable importance
    imp <- do.call(rbind,
                   lapply(names(mods), function(a){
                     
                     x <- mods[[a]]@variable.importance
                     x <- x/sum(x)
                     
                     data.frame(
                       studyID = studyID,
                       ES = es,
                       Algorithm = a,
                       Variable = names(x),
                       Importance = as.numeric(x)
                     )
                     
                   })
    )
    imp_all[[as.character(es)]] <- imp
    
    
    
  }
  
  imp_all <- dplyr::bind_rows(imp_all)
  
  performance_all <- dplyr::bind_rows(performance_all)
  row.names(performance_all)<-NULL
  
  saveRDS(
    imp_all,
    file.path(
      out_dir,
      paste0(studyID, "_importance.rds")
    )
  )
  
  saveRDS(
    performance_all,
    file.path(
      out_dir,
      paste0(studyID, "_algorithm_performance.rds")
    )
  )
  
  return(performance_all)
  
}



studID<-c("TRD")
#studID<-c("TRD")


#load stud areas
tmpSite<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",studID,"/raw_data_backup/study_site.gpkg"))

A_roi<-as.numeric(st_area(tmpSite))
resolution = 250^2
all_back_pts<- round(A_roi/resolution,0)
ind_pols<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",studID,"/raw_data_backup/all_ind_pols.gpkg"))%>%filter(siteID == studID)
## join imp_access to ind_pols via esID and userID from esmappingR1.csv
es_mapping<-read.csv(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",studID,"/raw_data_backup/es_mappingR1.csv"))
ind_pols <- ind_pols %>%
  dplyr::left_join(es_mapping %>% dplyr::select(esID, userID, imp_acc, confidence),
                   by = c("esID", "userID"))


esids<-unique(ind_pols$esID)

## site specific predictors
pred<-SSDM::load_var(path=paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",studID,"/raw_data_backup/2_env_var"), categorical = "lulc")


algo_comp<-run_algorithm_comparison(es_ids = esids,
                                    studyID = "TRD",
                                    ind_pols = ind_pols,
                                    pred = pred,
                                    A_roi = A_roi,
                                    all_back_pts = all_back_pts,
                                    out_dir = paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",studID,"/raw_data_backup"))
