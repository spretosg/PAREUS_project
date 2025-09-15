library(sf)
library(dplyr)

res<-st_read("C:/Users/reto.spielhofer/Downloads/SK021_PCA_landscape.geojson")
main_dir<-"P:/312204_pareus/"
sid<-"SK021"
stud_area<-read_sf(paste0(main_dir,"WP2/PGIS_ES_mapping/raw_data_backup/",sid,"/study_site.gpkg"))%>%filter(siteID %in% sid)
target_crs <- st_crs(stud_area)$wkt

#pa old
pa_strict<-res%>%filter(iucn_cat =="Ia" | iucn_cat == "II")
plot(pa_strict)


other_iucn<-res%>%filter(iucn_cat =="III" | iucn_cat == "IV"| iucn_cat == "V")

# new pa
pa_new<-res%>%filter(!is.na(rel_imp_pa))
library(ggplot2)


## centroids of new pa
path_new<-res%>%filter(least_cost_path_newPA > 1)

# Ensure same CRS
crs_target <- st_crs(stud_area)
pa_strict <- st_transform(pa_strict, crs_target)
pa_new    <- st_transform(pa_new, crs_target)

ggplot() +
  geom_sf(data = pa_strict, aes(fill = "Existing Protected Areas"), color = NA, alpha = 0.7) +
  # geom_sf(data = pa_new,    aes(fill = "New Protected Areas"),      color = NA, alpha = 0.7) +
  # geom_sf(data = path_new,    aes(fill = "OECM new"),      color = NA, alpha = 0.3) +
  # geom_sf(data = other_iucn,    aes(fill = "OECM existing"),      color = NA, alpha = 0.3) +
  geom_sf(data = stud_area, fill = NA, color = "black", size = 0.7) +
  scale_fill_manual(
    name   = "Protected Areas",
    values = c("Existing Protected Areas" = "#2ca02c")
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", axis.title = element_blank())
