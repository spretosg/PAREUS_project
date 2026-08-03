required <- c(
  "code",
  "data",
  "docs",
  "outputs",
  "README.md"
)

missing <- required[!file.exists(required)]

if(length(missing) > 0){
  
  stop(
    "Project root not found.\nMissing:\n",
    paste(missing, collapse="\n")
  )
  
}

message("✓ Project structure OK")
