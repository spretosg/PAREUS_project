library(readr)

catalogues <- list.files(
  "data",
  pattern = "datasets.csv",
  recursive = TRUE,
  full.names = TRUE
)

for(cat in catalogues){
  
  message("Reading ", cat)
  
  db <- read_csv(cat, show_col_types = FALSE)
  
  folder <- dirname(cat)
  
  for(i in seq_len(nrow(db))){
    
    outfile <- file.path(folder, db$filename[i])
    
    if(file.exists(outfile)){
      
      message("  ✓ ", db$filename[i])
      
      next
      
    }
    
    download.file(
      db$url[i],
      outfile,
      mode = "wb"
    )
    
  }
  
}