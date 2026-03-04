## script to calc stats of ES mapping

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

if(eval_round == "R1"){
  target_dir<-"/3_ind_R1"
}else{
  target_dir <-"/5_ind_R2"
}




# List the subfolders
subfolders <- list.dirs(paste0(main_dir,target_dir), full.names = TRUE, recursive = FALSE)

#import the mapping stats
map_stats<-read.csv(paste0(main_dir,"/es_mappingR1.csv"))
confidence<-map_stats%>%filter(siteID == stud_id )%>%select(confidence, esID)%>%group_by(esID)%>%summarise(med_conf =median(confidence))



### ---- in depth single es analysis ----
es<-"recr"
tar_folder<-subfolders[grepl(es, subfolders)]

raster_files <- list.files(tar_folder, full.names = TRUE)

# 2. Read all rasters
rasters <- lapply(raster_files, rast)

# 3. Assign names based on file names (for map layer labels)
layer_names <- tools::file_path_sans_ext(basename(raster_files))

# Step 3: Convert each raster to a data frame and add layer name
raster_dfs <- lapply(seq_along(rasters), function(i) {
  df <- as.data.frame(rasters[[i]], xy = TRUE)
  names(df)[3] <- "ES_capacity"  # Ensure the value column is named consistently
  df$layer_name <- layer_names[i]
  df
})

# Step 4: Combine all into one data frame
raster_all <- bind_rows(raster_dfs)

# Step 5: Plot using ggplot2
ggplot(raster_all, aes(x = x, y = y, fill = ES_capacity)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_equal() +
  facet_wrap(~ layer_name) +
  theme_minimal() +
  labs(fill = "ES capacity", title = "")

# variable importance

