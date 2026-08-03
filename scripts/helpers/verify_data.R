library(tools)

catalogues <- list.files(
  "data",
  pattern="datasets.csv",
  recursive=TRUE,
  full.names=TRUE
)

for(cat in catalogues){
  
  db <- read.csv(cat)
  
  folder <- dirname(cat)
  
  for(i in seq_len(nrow(db))){
    
    file <- file.path(folder, db$filename[i])
    
    if(!file.exists(file))
      stop(file, " missing.")
    
    if(nzchar(db$md5[i])){
      
      md5 <- md5sum(file)
      
      if(md5 != db$md5[i])
        warning(file, " checksum differs.")
      
    }
    
  }
  
}

message("✓ Data verified")