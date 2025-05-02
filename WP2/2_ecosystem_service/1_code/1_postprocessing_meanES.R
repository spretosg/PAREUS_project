### calculate mean raster per ES

library(terra)

main_dir<-"P:/312204_pareus/WP2/PGIS_ES_mapping/raw_data_backup/SK021"
eval_round<-"R1" #R2

if(eval_round == "R1"){
  target_dir<-"/3_ind_R1"
}else{
  target_dir <-"/5_ind_R2"
}

library(terra)


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

  # Create output filename
  folder_name <- basename(folder)
  out_file <- file.path(folder, paste0(folder_name, "_mean.tif"))

  # Save the mean raster
  writeRaster(mean_raster, out_file, overwrite = TRUE)

  cat("Saved mean raster to:", out_file, "\n")
}
