## optimization workflow with "prioritizr"
library(prioritizr)
library(terra)
library(dplyr)
#install.packages("highs", repos = "https://cran.rstudio.com/")

siteID<-"SK021"
main_dir<-"P:/312204_pareus/"

## import pu with costs
pu<-terra::rast(paste0(main_dir,"WP4/cost_raster_es/",siteID,"_optim_cost.tif"))
## import features
feat<-terra::rast(paste0(main_dir,"WP4/features/",siteID,"_optim_features.tif"))
print(feat)
# plot the first nine features
plot(feat[[1:2]], nr = 2, axes = FALSE)

## existing PA's for potential lock in
locked_in<-terra::rast(paste0(main_dir,"WP4/pa_existing/",siteID,"_pa.tif"))


## budget and problem
# calculate budget
budget <- terra::global(pu, "sum", na.rm = TRUE)[[1]] * 0.05

# create problem
p0 <-
  problem(pu, features = feat) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.05) %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0.1, verbose = FALSE)

# print problem
print(p0)

# solve the problem
# lock out the already protected areas
p1 <-
  p0 %>%
  add_locked_in_constraints(locked_in)

# solve the problem
s1 <- solve(p1)

# plot the solution
#plot(s2, main = "Solution", axes = FALSE)

# plot solutions
plot(
  c(s0, s1), main = c("baseline", "no PA"),
  axes = FALSE, type = "classes", col = c("grey70", "darkgreen")
)

# extract the objective
#print(attr(s1, "objective"))
# extract state message from the solver
print(attr(s1, "status"))

# plot the solution
plot(s0, main = "Solution 0", axes = FALSE)



### evaluation
# calculate number of selected planning units by solution
eval_n_summary(p1, s1)

# calculate total cost of solution
eval_cost_summary(p1, s1)


# calculate target coverage for the solution
p1_target_coverage <- eval_target_coverage_summary(p1, s1)
print(p1_target_coverage)

# check percentage of the features that have their target met given the solution
print(mean(p1_target_coverage$met) * 100)


### new version with existing pa's
