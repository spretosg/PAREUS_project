### downloading and storage from gcs to local of all data
library(bigrquery)
library(dplyr)
library(sf)
# library(SSDM)
# library(tidyverse)
library(googleCloudStorageR)
# library(rnaturalearth)
# library(rnaturalearthdata)

#store the raw data on nina servers
stud_id<-"FRA_BAR2"
out_master_path<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",stud_id,"/raw_data_backup")
dev<-"dev"
id<-"pareus"

#############


dataset<-paste0(id,"_",dev)

if(id == "wendy"){
  project_id<-"eu-wendy"
}else{
  project_id<-id
}

bq_auth(
  path = paste0("docs/",project_id,"_key.json")
)


cred_path<-paste0("docs/",project_id,"_key.json")
gcs_auth(cred_path)
bucket<-paste0(id,"_geopros_",dev)
gcs_global_bucket(bucket)

con_admin<-data.frame(
  project = project_id,
  dataset = dataset,
  billing =project_id
)


con_admin <- dbConnect(
  bigrquery::bigquery(),
  project = con_admin$project,
  dataset = con_admin$dataset,
  billing = con_admin$billing
)

tables <- bq_dataset_tables(bq_dataset(project_id, dataset))
tables<-tables[-1]

### store all bq tables as csv or gpkg if not present
if(is.null(list.files(out_master_path))){
  for (table_ref in tables) {
    table_id <- table_ref$table
    full_table <- bq_table(project_id, dataset, table_id)

    # Get schema
    schema <- bq_table_meta(full_table)$schema$fields
    field_names <- sapply(schema, function(f) f$name)
    field_types <- sapply(schema, function(f) f$type)

    # Check if there's a GEOGRAPHY column
    has_geo <- any(field_names == "geometry")

    # Read the table
    df <- bq_table_download(full_table)

    if (has_geo) {
      # Convert GEOGRAPHY column to sf geometry
      geo_col <- field_names[field_names == "geometry"][1]
      df_sf <- st_as_sf(df, wkt = geo_col, crs = 4326)
      out_path <- file.path(out_master_path, paste0(table_id, ".gpkg"))
      st_write(df_sf, out_path, delete_dsn = TRUE, quiet = TRUE)
      message(paste("Saved", table_id, "as GPKG"))
    } else {
      out_path <- file.path(out_master_path, paste0(table_id, ".csv"))
      write.csv(df, out_path)
      message(paste("Saved", table_id, "as CSV"))
    }
  }
}

## download all raster:
pattern<-paste0(stud_id,"/3_ind_R1")

objects <- gcs_list_objects(prefix = pattern)
# Create local folder to store the files
local_dir <- paste0(out_master_path,"/",pattern)
dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)

# Loop over and download each object
for (obj in objects$name) {
  message("Downloading: ", obj)

  # Strip the prefix from path to preserve folder structure locally
  relative_path <- sub(paste0("^",stud_id,"/5_ind_R2"), "", obj)
  local_path <- file.path(local_dir, relative_path)
  dir.create(dirname(local_path), recursive = TRUE, showWarnings = FALSE)

  gcs_get_object(obj, saveToDisk = local_path, overwrite = TRUE)
}

