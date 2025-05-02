## prep cost layers based on ecosystem services according to cookbook Nr.4
library(terra)
library(sf)
library(dplyr)
library(ppcor)


main_dir<-"P:/312204_pareus/WP2/PGIS_ES_mapping/raw_data_backup/SK021"
## Read mean es raster from R2 of workshop
erosion = rast(list.files(paste0(main_dir,"/3_ind_R1/erosion"), full.names = TRUE))
mat  = rast(list.files(paste0(main_dir,"/3_ind_R1/mat"), full.names = TRUE))
farm  = rast(list.files(paste0(main_dir,"/3_ind_R1/farm"), full.names = TRUE))
flood  = rast(list.files(paste0(main_dir,"/3_ind_R1/flood"), full.names = TRUE))
wild_plant  = rast(list.files(paste0(main_dir,"/3_ind_R1/wild_plant"), full.names = TRUE))
wild_hunt  = rast(list.files(paste0(main_dir,"/3_ind_R1/wild_hunt"), full.names = TRUE))
sense  = rast(list.files(paste0(main_dir,"/3_ind_R1/sense"), full.names = TRUE))


es<-c("erosion","mat","farm","flood","wild_plant","wild_hunt","sense")
# from AHP in tool
w<-c(0.2,0.4,0.1,0.3)
es<-data.frame(es,w)
colnames(es)<-c("es","weight")

## define the pcor function using three rasters of spec dimensions --> here 20 individual ES1 rasters
# n_part<-nlyr(aest)
# fun_pcor =  function(x) {
#   # n_part=20
#   Rs = pcor.test(x[1:n_part],x[(n_part+1):(2*n_part)],x[(2*n_part+1):(3*n_part)])
#   Rx = Rs$estimate
#   # Px = Rs$p.value
#   # return(c(Rx, Px))
#   return(Rx)
# }

n_es<-as.numeric(nrow(es))

fun_pcor_alt =  function(x) {
  # n_part=20
  #first two are the main, the others for the correction, take the total number of es to select the correct array values
  Rs = pcor.test(x[[1:4]],x[(n_part+1):(2*n_part)],x[(2*n_part+1):(n_es*n_part)],method = "spearman")
  Rx = Rs$estimate
  # Px = Rs$p.value
  # return(c(Rx, Px))
  return(Rx)
}

## generate all combinations but not ES1 == ES1...

comb<-subset(expand.grid(rep(list(es$es),2)), Var1 != Var2)

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
  nlayers <- min(nlyr(get(left)), nlyr(get(right)))
  rank1 <- app(get(left), rank)
  rank2 <- app(get(right), rank)
  right <- rank1[[1:nlayers]]
  left <- rank2[[1:nlayers]]
  
  sji <- app(c(rank1, rank2), fun = function(x) {
    n <- length(x) / 2
    cor(x[1:n], x[(n+1):(2*n)], method = "pearson")
  })
  ## exclusiveness of an es is then the 1-Sji raster (values close to 1 == exclusive) however pos or neg correlation does not matter
  eji<-1-abs(sji)
  ## multiply exclusiveness raster with benefit (mean of all participants) raster i
  bi<-mean(left)
  bj<-mean(right)
  wj<-as.numeric(es%>%filter(es==as.character(a$Var2))%>%dplyr::select(weight))
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
for(n in 1: nrow(es)){
  tmp_es<-es[n,]$es
  pattern <- paste0(tmp_es,"_")
  selected_rasters <- c_stack[[grep(pattern, names(c_stack))]]
  ci<-sum(selected_rasters)*as.numeric(es%>%filter(es==tmp_es)%>%dplyr::select(weight))
  cost<-c(cost,ci)
}
#summarize the cost of individual es i
c=sum(cost)
plot(c)

#write raster
terra::writeRaster(c,paste0(main_dir,"/cost_raster_es.tif"))




################



# Ensure equal number of layers: match by time, interpolate, or truncate


plot(eji, main = "Spearman Rank Correlation (per pixel)")
