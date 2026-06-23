library(circuitscaper)
#cs_install_julia()
Sys.which("julia")
system("julia --version")
#cs_setup(threads = 4)

# library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
# library(gdistance)

source("WP4/1_code/wp4_functions_utils.R")

main_dir<-"P:/312204_pareus/"
siteID<-"FRL04"

stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/study_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="FRL04")

pu<-st_read(paste0("WP4/2_output/02_optim/",siteID,"_input_final_grid.json"))# pu<-st_read("WP4/2_output/02_optim/FRL04_input_final_grid.geojson")


es_cond<-paste0(main_dir,"WP4/features/",siteID,"_ec.tif")
es_cond<-rast(es_cond)
es_cond

es_cond<-project(es_cond,st_crs(pu)$wkt) 

resistance<-log(1/es_cond )
#resistance <- 1/es_cond
resistance
### connectivity

# calculate connectivity of optimized core protected areas
# cells <- pu %>% 
#   filter(solution_1>0 | lock_in == T)
# 
# # 2. Union touching cells
# core_pa <- cells %>%
#   summarise(geometry = st_union(geometry)) %>%
#   st_cast("POLYGON")%>%st_transform(st_crs(resistance))
# core_cent <- st_coordinates(st_centroid(core_pa))
# 
r_coarse <- aggregate(resistance, fact = 2)
r_coarse

# 
# pw_result <- cs_pairwise(r_coarse, core_cent)
# plot(pw_result$current_map)
# plot(st_centroid(core_pa),add=T, col= "red")
# writeRaster(pw_result$current_map,paste0("WP4/2_output/02_optim/pw_connectivity_",siteID,".tif"))
r_coarse[r_coarse <= 0] <- 1e-6
plot(r_coarse)

# calculate moving window connectivity based on ecosystem condition in landscape
start<-Sys.time()
mw_result <- os_run(r_coarse, radius = 20, block_size = 10)
plot(mw_result$normalized_current)
print(Sys.time()-start)
# crs(result2) <- as.character(crs(r_coarse))
mw_result <- terra::project(mw_result, sf::st_crs(pu)$wkt)
writeRaster(mw_result,paste0("WP4/2_output/02_optim/mw_connectivity_",siteID,".tif"))

# mw_result<-rast(paste0("WP4/2_output/02_optim/mw_connectivity_",siteID,".tif"))

# sample connetivity values to pu
pu$mw_connectivity<- terra::extract(mw_result$normalized_current, pu, fun = mean, na.rm = TRUE)[,2]

conn<-ggplot(pu) +
    geom_sf(aes(fill = mw_connectivity), color = NA) +
    scale_fill_viridis_c(option = "plasma", name = "Connectivity") +
    geom_sf(data = stud_area, fill = NA, color = "black") +
    theme_minimal()+
    theme(text = element_text(size = 20))
  #ggsave(paste0("WP4/2_output/02_optim/",siteID,"_OECM_suit.png"), plot = p, width = 8, height = 6, dpi = 300)


pu<- zero_one_scale(
  pu,
  cols = c("mw_connectivity")
)

pu<-pu%>%mutate(LULC_class = case_when(sampled_habitat ==1 ~ "artificial_surfaces",
                                       sampled_habitat ==2 ~ "agricultural_areas",
                                       sampled_habitat ==3 ~ "forests",
                                       sampled_habitat ==4 ~ "wetlands",
                                       sampled_habitat ==5 ~ "water"))

ggplot(pu, aes(x = LULC_class, y = mw_connectivity, fill = LULC_class, group = LULC_class)) +
  geom_boxplot() +
  # scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  theme_minimal() +
  theme(legend.position = "none",text = element_text(size = 20))+
  labs(
    x = "",
    y = "Connectivity",
    fill = ""
  )



