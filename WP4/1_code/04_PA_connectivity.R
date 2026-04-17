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
siteID<-"SK021"

stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/stud_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="SK021")

pu<-st_read(paste0("WP4/2_output/02_optim/PA_optim_",siteID,".geojson"))
# pu<-st_read("WP4/2_output/02_optim/FRL04_input_final_grid.geojson")


es_cond<-paste0(main_dir,"WP4/features/",siteID,"_ec.tif")
es_cond<-rast(es_cond)
es_cond

es_cond<-project(es_cond,st_crs(pu)$wkt) 

resistance<-log(1/es_cond )
resistance <- 1/es_cond
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
crs(result2) <- as.character(crs(r_coarse))
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



### OECM suitability
w_cult = 0.25
w_prov = 0.25
w_ec = 0.25
w_connect = 0.25

pu$oecm_suit<-oecm_lin_w(pu$sampled_cult_scaled,pu$sampled_prov_scaled,
                         pu$sampled_condition_scaled,pu$mw_connectivity_scaled,w_cult,w_prov,w_ec,w_connect)

suit<-ggplot(pu) +
  geom_sf(aes(fill = oecm_suit), color = NA) +
  scale_fill_viridis_c(option = "viridis", name = "OECM suitability") +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 20))

ggplot(pu, aes(x = LULC_class, y = oecm_suit, fill = LULC_class, group = LULC_class)) +
  geom_boxplot() +
  # scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  theme_minimal() +
  theme(legend.position = "none",text = element_text(size = 20))+
  labs(
    x = "",
    y = "OECM suitability",
    fill = ""
  )

maps<-ggarrange(conn, suit)
ggsave(paste0("WP4/2_output/02_optim/",siteID,"_OECM_connect.png"), plot = maps, width = 18, height = 10, dpi = 300)


p_oecm<-0.2 # from total area

#exclude existing and new core pa here
oecm <- pu %>% filter(oecm_suit >= quantile(oecm_suit, 1-p_oecm, na.rm = TRUE))

oecm <- oecm %>% filter(is_core_prot_lulc == F)
core_pa <- pu %>% filter(is_core_prot_lulc == T)

other_pa<-pu%>%filter(max_IUCN_class<6)


## highest 20% OECM suitability

oecm$pa_group <- ifelse(oecm$n_pa > 0, "High OECM suitability already protected (IUCN III-VI)", "High OECM suitability not protected - potential OECM")

## other pa, remove ids from OECM
other_pa_filt <- other_pa[!other_pa$ID %in% oecm$ID, ]

other_pa_filt <- other_pa_filt[!other_pa_filt$ID %in% core_pa$ID, ]
other<-pu %>% filter(core_pa_lulc == "other")
other <- other[!other$ID %in% oecm$ID, ]
other <- other[!other$ID %in% other_pa_filt$ID, ]

stats_oecm<-oecm%>%group_by(pa_group,sampled_habitat)%>%summarise(area = sum(area)/10^6)%>%st_drop_geometry()
stats_other_pa<-other_pa_filt%>%group_by(sampled_habitat)%>%summarise(area = sum(area)/10^6)%>%st_drop_geometry()
stats_core_pa<-core_pa%>%group_by(sampled_habitat)%>%summarise(area = sum(area)/10^6)%>%st_drop_geometry()
stats_other<-other%>%group_by(sampled_habitat)%>%summarise(area = sum(area)/10^6)%>%st_drop_geometry()

library(ggnewscale)
cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#00FF00",
  "proposed upgrade existing PA" = "#B8FFB8"
)

plot(st_geometry(oecm),col="purple",add = T)
plot(st_geometry(core_pa),col="red",add=T)

plot(st_geometry(other_pa_filt),col="lightblue",add=T)
plot(st_geometry(other),col="lightgrey")

p<-ggplot() +
  # Layer 0: other pa
  geom_sf(data = other_pa,
          aes(fill = "Other protected areas (IUCN III-VI)"),
          color = NA) +
  
  # Layer 1: oecm (purple based on n_pa)
  geom_sf(data = oecm,
          aes(fill = pa_group),
          color = NA) +
  scale_fill_manual(
    values = c("High OECM suitability already protected (IUCN III-VI)" = "#E6D9FF",
               "High OECM suitability not protected - potential OECM" = "#4B0082",
               "Other protected areas (IUCN III-VI)" = "#ADD8E6"),
    name = NULL
  )  +
  
  new_scale_fill() +
  
  # Layer 2: core_pa (categorical colors)
  geom_sf(data = core_pa,
          aes(fill = core_pa_lulc),
          color = NA) +
  scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  new_scale_fill() +
  
  # stud area
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()
ggsave(paste0("WP4/2_output/02_optim/",siteID,"_optim_PCA_landscape.png"), plot = p, width = 12, height = 16, dpi = 300)
