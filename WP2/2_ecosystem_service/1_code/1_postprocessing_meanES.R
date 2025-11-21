### calculate mean raster per ES

library(terra)
source("WP2/wp2_functions_utils.R")
stud_id<-"FRA_BAR2"
main_dir<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",stud_id,"/raw_data_backup")
eval_round<-"R1" #R2

if(eval_round == "R1"){
  target_dir<-"/3_ind_R1"
}else{
  target_dir <-"/5_ind_R2"
}




# List the first 10 subfolders
subfolders <- list.dirs(paste0(main_dir,target_dir), full.names = TRUE, recursive = FALSE)[1:10]

# Loop through each subfolder
for (folder in subfolders) {

  # List all raster files (you can adjust the pattern if needed)
  raster_files <- list.files(folder, full.names = TRUE)

  # Read all rasters into a SpatRaster list
  rasters <- lapply(raster_files, rast)

  # Stack the rasters (assuming same extent/resolution)
  raster_stack <- rast(rasters)

  # Calculate mean raster
  mean_raster <- mean(raster_stack, na.rm = TRUE)
  cv_raster <- cv_rast(raster_stack)

  # Create output filename
  # folder_name <- basename(folder)
  mean_path<-paste0(main_dir,"/4_mean_",eval_round)
  cv_path<-paste0(main_dir,"/5_cv_",eval_round)
  if(!exists(mean_path)){
    dir.create(mean_path)
  }
  if(!exists(cv_path)){
    dir.create(cv_path)
  }
  out_mean<- file.path(mean_path, paste0(basename(folder), "_mean.tif"))
  out_cv<- file.path(cv_path, paste0(basename(folder), "_cv.tif"))
  names(mean_raster)<-basename(folder)
  names(cv_raster)<-basename(folder)

  writeRaster(mean_raster, out_mean, overwrite = TRUE)
  writeRaster(cv_raster, out_cv, overwrite = TRUE)

  cat("Saved to:", main_dir, "\n")
}
