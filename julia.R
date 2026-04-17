# 1. Install the R package
remotes::install_github("matthewkling/circuitscaper")

# 2. Let the package install Julia and the necessary Julia libraries
library(circuitscaper)
cs_install_julia()
