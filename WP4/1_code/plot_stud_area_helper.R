
library(sf)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(ggrepel)

main_dir<-"P:/312204_pareus/"
siteID<-"FRL04"
stud_areas<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/study_site.gpkg"))
targetSites<-c("SK021","FRL04","FRA_BAR2")
nor<-c("NOR_lowlands","NOR_coast","NOR_mtn")

others<-stud_areas%>%filter(siteID %in% targetSites)%>%select(cntrID,siteID,siteNAME)
nor_areas<-stud_areas%>%filter(siteID %in% nor)%>%
  group_by(cntrID)%>%summarise()
nor_areas$siteID = "NOR_TRD"
nor_areas$siteNAME = "Trøndelag"
stud_areas<-rbind(others,nor_areas)

europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

stud_areas <- st_transform(stud_areas, st_crs(europe))

coords <- st_coordinates(st_centroid(stud_areas)) %>%
  as.data.frame() %>%
  cbind(st_drop_geometry(stud_areas))


ggplot() +
  geom_sf(data = europe, fill = "gray95", color = "gray50") +
  geom_sf(data = stud_areas, fill = "red", alpha = 0.6) +
  geom_text_repel(
    data = coords,
    aes(X, Y, label = siteID),
    nudge_x = c(0.5,-1,0.3,2),
    nudge_y = c(1.1,-2,1.5,2.5)
  ) +
  coord_sf(xlim = c(-15, 40), ylim = c(35, 75), expand = FALSE) +
  theme_minimal() 
