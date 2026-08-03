#wp2 functions

cv_rast<-function(r){
  cv=terra::stdev(r,na.rm=T)/mean(r,na.rm=T)
  return(cv)
}