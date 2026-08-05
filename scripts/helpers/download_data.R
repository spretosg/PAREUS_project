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

## create intermediate storage folders
dirs <- c(
  "outputs",
  "outputs/WP1",
  "outputs/WP2",
  "outputs/WP3",
  "outputs/WP4",
  "outputs/WP4/01_PA_analysis",
  "outputs/WP4/02_optim",
  "outputs/WP4/03_pca_landscape"
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}
