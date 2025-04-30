## prep cost layers based on ecosystem services according to cookbook Nr.4
library(terra)
library(sf)
library(dplyr)
library(ppcor)

# cntr<-"NO"
##1. Ecosystem service benefits (Bi)
### these are the maps resulting from es_accounting and/or geoprospective approach
# dat_path<-"C:/Users/reto.spielhofer/git/PAREUS_routines/2_optimal_siting_pa"
# es_files<-list.files(paste0(dat_path,"/2_protected_features"),full.names = T)
# 
# es_raster <- rast(es_files[1:3])



### test pcor
# create sample data frame

### Synergy between ES

## here we read the individual es rasters for e.g. N=20 participants
es1  = rast(array(runif(10*10*20),c(10,10,20)))
es2  = rast(array(runif(10*10*20),c(10,10,20)))
es3  = rast(array(runif(10*10*20),c(10,10,20)))
es4  = rast(array(runif(10*10*20),c(10,10,20)))
es5  = rast(array(runif(10*10*20),c(10,10,20)))


es<-c("es1","es2","es3","es4","es5")
# from AHP in tool
w<-c(0.2,0.4,0.1,0.2,0.1)
es<-data.frame(es,w)
colnames(es)<-c("es","weight")

## define the pcor function using three rasters of spec dimensions --> here 20 individual ES1 rasters
n_part<-nlyr(es1)
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
  Rs = pcor.test(x[1:n_part],x[(n_part+1):(2*n_part)],x[(2*n_part+1):(n_es*n_part)],method = "spearman")
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
  left<-as.character(a$Var1)
  right<-as.character(a$Var2)
  # select the others
  ctrl<-es%>%filter(es != left & es != right)%>%dplyr::select(es)
  # empty stack for the r raster of the es i left
  # r_stack <- rast()
  # 
  # for(j in 1:nrow(no_comb)){
  #   z = c(get(left),get(right),get(as.character(no_comb[j,])))
  #   
  #   # z = c(es1,es2,es3)
  #   pcor_temp<-app(z, fun_pcor)
  #   # write in stack
  #   # p_stack<-raster::stack(p_stack,pcor_temp$Px)
  #   r_stack<-c(r_stack,pcor_temp)
  # }
  
  ## alternative make a vector containing all the raster in the correct order for the pcor_alt function:
  
  base<-c(get(left),get(right))
  control_rast<-c()
  for(n in 1: nrow(ctrl)){
    control_rast<-c(control_rast,get(as.character(ctrl[n,])))
  }
  control_rast<-rast(control_rast)
  z<-c(base,control_rast)
  
  # z = c(get(left),get(right),get(as.character(no_comb[1,])),get(as.character(no_comb[2,])))
  sji<-app(z, fun_pcor_alt)
  
  #the mean of all r rasters represents the synergy-trade off between Esi and ESj given all other ES
  # sji<-mean(r_stack,na.rm=T)
  # rm(r_stack)
  
  ## exclusiveness of an es is then the 1-Sji raster (values close to 1 == exclusive) however pos or neg correlation does not matter
  eji<-1-abs(sji)
  ## multiply exclusiveness raster with benefit (mean of all participants) raster i
  bi<-mean(get(left))
  bj<-mean(get(right))
  wj<-as.numeric(es%>%filter(es==right)%>%dplyr::select(weight))
  c_rast<-eji*bi*bj*wj
  names(c_rast)<-paste0(left,"_",right)
  
  names(sji)<-paste0(left,"_",right)
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







