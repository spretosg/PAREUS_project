## prep cost layers based on ecosystem services according to cookbook Nr.4
library(terra)
library(sf)
library(dplyr)
library(ppcor)

stud_id<-"FRL04"
main_dir<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",stud_id,"/raw_data_backup")

out_dir<-"P:/312204_pareus/WP4/cost_raster_es"
## Read mean es raster from R2 of workshop

mean_rast<-list.files(paste0(main_dir,"/4_mean_R1"), full.names = TRUE)
mean_rast<-rast(mean_rast)

es<-names(mean_rast)
# from AHP in tool
ahp_w<-read.csv(paste0(main_dir,"/ahp_weights.csv"))


n_es<-as.numeric(length(es))


comb<-subset(expand.grid(rep(list(es),2)), Var1 != Var2)

# empty stack for synergy s and cost c
s_stack<-rast() # stack to store the synergies  between esi and esj controlled for all others
c_stack <- rast() # stack to store the costs of esi

## for the length in comb we need to choose the two es and the other es to pcor with
for(i in 1: nrow(comb)){
  a<-comb[i,]
  print(i/nrow(comb))
  left<-as.character(a$Var1)
  right<-as.character(a$Var2)
  # select the others
  #nlayers <- min(nlyr(get(left)), nlyr(get(right)))
  # rank1 <- app(get(left), rank)
  # rank2 <- app(get(right), rank)
  # right <- rank1[[1:nlayers]]
  # left <- rank2[[1:nlayers]]
  # 
  # sji <- app(c(right, left), fun = function(x) {
  #   n <- length(x) / 2
  #   #print(n)
  #   cor(x[1:n], x[(n+1):(2*n)], method = "pearson")
  # })
  l_raster <- mean_rast[[left]]
  r_raster <- mean_rast[[right]]
  tmp_cor_rast<-c(l_raster,r_raster)
  ## local, focal correlation with 5x5 window
  sji<-focalPairs(tmp_cor_rast,matrix(1,5,5),"pearson")
  
  ## exclusiveness of an es is then the 1-Sji raster (values close to 1 == exclusive) however pos or neg correlation does not matter
  eji<-1-abs(sji)
  ## multiply exclusiveness raster with benefit (mean of all participants) raster i
  bi<-mean(values(l_raster),na.rm=T)
  bj<-mean(values(r_raster),na.rm=T)
  wj<-as.numeric(ahp_w%>%filter(esID==as.character(a$Var2))%>%dplyr::select(pref_adj))
  c_rast<-eji*bi*bj*wj
  names(c_rast)<-paste0(as.character(a$Var1),"_",as.character(a$Var2))
  
  names(sji)<-paste0(as.character(a$Var1),"_",as.character(a$Var2))
  ## assign es comb to both
  s_stack<-c(s_stack,sji)
  c_stack<-c(c_stack,c_rast)
}

# finally the global ES cost is given as the sum of all raster individual ES costs ci in the c_stack with left = i and multiply with the weight i
# Define your pattern
cost<-rast()
for(n in 1: length(es)){
  tmp_es<-es[n]
  pattern <- paste0(tmp_es,"_")
  selected_rasters <- c_stack[[grep(pattern, names(c_stack))]]
  ci<-sum(selected_rasters,na.rm=T)*as.numeric(ahp_w%>%filter(esID==tmp_es)%>%dplyr::select(pref_adj))
  cost<-c(cost,ci)
}
#summarize the cost of individual es i
c=sum(cost,na.rm=T)
plot(c)

#write raster
terra::writeRaster(c,paste0(out_dir,"/",stud_id,"_cost_raster_es.tif"),overwrite=T)

