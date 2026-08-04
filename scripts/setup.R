############################################################
## PAREUSD project setup
## Installs packages, downloads data and checks project
############################################################

message("=== Project setup ===")

## Step 1 --------------------------------------------------
source("scripts/helpers/check_project.R")

## Step 2 --------------------------------------------------
if (file.exists("renv.lock")) {
  
  message("Restoring R environment...")
  
  if (!requireNamespace("renv", quietly = TRUE))
    install.packages("renv")
  
  renv::restore(prompt = FALSE)
  
} else {
  
  message("No renv.lock found.")
}

## Step 3 --------------------------------------------------
source("scripts/helpers/download_data.R")

## Step 4 --------------------------------------------------
source("scripts/helpers/verify_data.R")

# ## Step 5 --------------------------------------------------
# source("scripts/helpers/download_optional.R")
# 
# ## Step 6 --------------------------------------------------
# source("scripts/helpers/check_project.R")

## Step 7 --------------------------------------------------
message("")
message("==============================")
message("Project successfully installed")
message("==============================")