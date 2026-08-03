library(sf)
library(leaflet)
library(DT)
library(dplyr)
library(ggplot2)
library(dplyr)
library(terra)

source("WP2/wp2_functions_utils.R")
main_dir<-"P:/312204_pareus/"


sites<-st_read(paste0(main_dir,"pareus_repository/pareus_sites.gpkg"))

sites<-st_make_valid(sites)

####----user statistics----
user_stats<-read.csv(paste0(main_dir,"pareus_repository/WP2/ES_mapping/mapper.csv"))


user_poly<-st_read(paste0(main_dir,"pareus_repository/WP2/ES_mapping/user_polygons.gpkg"))
user_poly<-st_make_valid(user_poly)
df1<-user_poly%>%st_drop_geometry()%>%group_by(siteID,esID)%>%summarise(n_part = n_distinct(userID))
ggplot(df1, aes(x = esID, y = n_part)) +
  geom_col(fill = "#4C78A8") +
  facet_wrap(~siteID, ncol = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "Ecosystem service", y = "number of participants")

df2<-user_poly%>%group_by(siteID,esID)%>%summarise(sum_area = as.numeric(sum(st_area(geom))/10^6),
                                                   average_part_area = sum_area/n_distinct(userID))%>%st_drop_geometry()

ggplot(df2, aes(x = esID, y = average_part_area)) +
  geom_col(fill = "grey") +
  facet_wrap(~siteID, ncol = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "Ecosystem service", y = "Mapped area / participant [km2]")

####---- importance of ES ----
es_imp<-read.csv(paste0(main_dir,"pareus_repository/WP2/ES_mapping/weights_all.csv"),sep=";")
ggplot(es_imp, aes(x = esID, y = pref_adj)) +
  geom_col(fill = "grey") +
  facet_wrap(~siteID, ncol = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "Mean indicated importance")
####---- single ES ----


####---- mean es maps ----
mean_files <- list.files(
  paste0(main_dir,"/pareus_repository/WP2/ES_mapping/"),
  pattern = "_Wmean\\.tif$",
  full.names = TRUE
)
studyIDs <- sub(
  "^(.*)_Wmean$",
  "\\1",
  tools::file_path_sans_ext(basename(mean_files))
)

plots <- lapply(studyIDs, function(st){
  
  f <- mean_files[grepl(paste0( st,"_Wmean"), basename(mean_files))]
  tmp_area<-sites%>%filter(siteID == st)
  
  r <- rast(f)
  r <- mask(crop(r, vect(tmp_area)), vect(tmp_area))
  
  df <- as.data.frame(r , xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "es_mean"
  
  df$studyID <- st
  
  ggplot(df) +
    geom_raster(aes(x, y, fill = es_mean)) +
    scale_fill_viridis_c("turbo") +
    geom_sf(
      data = tmp_area,
      fill = NA,
      colour = "black",
      linewidth = 0.3
    )+
    theme_bw() +
    ggtitle(st) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  
})


# Combine into one figure
combined_plots <- patchwork::wrap_plots(plots,
                                        ncol = 2) +
  patchwork::plot_layout(guides = 'collect')

combined_plots
