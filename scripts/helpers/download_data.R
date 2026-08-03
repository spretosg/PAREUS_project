library(readr)
library(dataverse)

Sys.setenv(DATAVERSE_SERVER = "dataverse.harvard.edu")

catalogue <- read.csv("data_catalogue.csv", sep=";")

for (k in seq_len(nrow(catalogue))) {
  
  message("Downloading ", catalogue$WP[k])
  
  ds <- get_dataset(catalogue$DOI[k])
  
  dir.create(catalogue$Destination[k],
             recursive = TRUE,
             showWarnings = FALSE)
  
  for (i in seq_len(nrow(ds$files))) {
    
    raw <- get_file(ds$files$id[i])
    
    writeBin(
      raw,
      file.path(
        catalogue$Destination[k],
        ds$files$label[i]
      )
    )
  }
}