# GAP analysis
#how much of each natural land cover type is protected through IUCN cat Ia and II?
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)
source("WP4/1_code/wp4_functions_utils.R")

# siteID<-"FRL04"
#siteID<-"SK021"
siteID<-"FRL04"
target_core_prot_fraction <-0.1 #how much of each biome should be protected strictly
main_dir<-"P:/312204_pareus/"

#inputs

## stud_area 
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/study_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="FRL04")
# target_crs <- 2154 #adjust this france
target_crs<-25833
stud_area<-st_transform(stud_area,target_crs)
total_area<-st_area(stud_area)
stud_ara_vect <- vect(stud_area)

grid <- st_make_grid(stud_area, cellsize = 1200, square = F) # depending on the size FRA SVK 750 TRD 1200
grid <- st_sf(geometry = grid)


grid<-st_intersection(grid,stud_area["siteID"])
grid$ID<-c(1:nrow(grid))
grid$area<-as.numeric(st_area(grid))

# extract habitat per pu 
lulc<-paste0(main_dir,"WP4/habitat/",siteID,"_lulc.tif")
lulc<-terra::rast(lulc)
lulc <- project(lulc, paste0("epsg:",target_crs))
lulc <- crop(lulc, stud_ara_vect)
lulc <- mask(lulc, stud_ara_vect)
lulc<-floor(lulc / 100)
lulc_highres <- disagg(lulc, fact = 5, method = "near")
lulc_highres <- crop(lulc_highres, stud_ara_vect)
lulc_highres <- mask(lulc_highres, stud_ara_vect)

grid$sampled_habitat <- terra::extract(
  lulc_highres,
  vect(st_centroid(grid))
)[,2]

#remove grid outside stud area
grid<-grid%>%filter(!is.na(sampled_habitat))

lulc_stats<-grid%>%group_by(sampled_habitat)%>%summarise(area_km2=sum(as.numeric(area))/10^6)

## sample PA status 
PA<-st_read(paste0(main_dir,"WP4/pa_existing/WDPA_",siteID,".shp"))
PA<-st_transform(PA,st_crs(target_crs))
PA <- st_intersection(PA, stud_area)
PA<-st_make_valid(PA)


PA <- PA %>%
  mutate(class = case_when(
    IUCN_CAT == "Not Applicable" ~ 1,
    IUCN_CAT == "Not Assigned" ~ 1,
    IUCN_CAT == "Not Reported" ~ 1,
    IUCN_CAT == "IV" ~ 2,
    IUCN_CAT == "V" ~ 3,
    IUCN_CAT == "IV" ~ 4,
    IUCN_CAT == "III" ~ 5,
    IUCN_CAT == "II" ~ 6,
    IUCN_CAT == "Ia" ~ 7,
    TRUE ~ 0
  ))

#crop to study area
PA<-st_intersection(PA,stud_area)

grid<-grid%>%mutate(lulc_name = case_when(sampled_habitat == 1 ~ "built-up",
                                          sampled_habitat == 2 ~ "agricultural",
                                          sampled_habitat == 3 ~ "forest",
                                          sampled_habitat == 4 ~ "wetland",
                                          sampled_habitat == 5 ~"water"))

#sampling number of PA in one cell
# get intersection index list
ints <- st_intersects(grid, PA)

# count unique IUCN categories per grid cell
grid$n_pa <- sapply(ints, function(i) {
  length(unique(PA$IUCN_CAT[i]))
})

#sample highest status of protection
grid <- st_join(grid, PA[, c("class", "IUCN_CAT")], join = st_intersects)
grid <- grid %>%
  group_by(ID) %>%
  slice_max(class, n = 1, with_ties = FALSE) %>%
  ungroup()

### core prot vs other pa's
grid<-grid%>%mutate(exisiting_corePA = case_when(class<6 ~ F, class>5 ~ T, is.na(class)~F))

gap_glob <- protection_gap(
  grid,
  lulc_col = "lulc_name",
  corePA_col = "exisiting_corePA",
  lockout = "built-up",
  mode = "global",
  target = 0.1
)

gap_lulc <- protection_gap(
  grid,
  lulc_col = "lulc_name",
  corePA_col = "exisiting_corePA",
  lockout = "built-up",
  mode = "class",
  target = c(forest = 0.1, water = 0.2, wetland = 0.2, agricultural = 0.05)
)

gap<-rbind(gap_lulc,gap_glob)
gap<-gap%>%filter(!is.na(lulc_name))

#overlaying PAs
# p<-ggplot() +
#      geom_sf(
#        data = grid%>%filter(n_pa>0),
#         aes(fill = n_pa),
#      color = "NA"  )+
#   scale_fill_gradient(
#     high = "#C50202",
#     low = "#FFC2C2",
#     space = "Lab",
#     na.value = "grey50",
#   )+
#   geom_sf(data = stud_area, fill = NA, color = "black") +
#   theme_minimal()+
#   theme(text = element_text(size = 20))
# ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_n_pa.png"), plot = p, width = 8, height = 6, dpi = 300)
# 

# PA classes
p<-ggplot() +
  geom_sf(
    data = grid,
    aes(fill = IUCN_CAT),
    color = "NA"  )+
  scale_fill_manual(
    values = c(
      "Ia" = "#00A300",
      "II" = "#00FF00",
      "V" = "#B589D6",
      "Not Applicable" = "#F1EE8D",
      "Not Assigned" = "#F1EE8D",
      "Not Reported" = "#F1EE8D",
      "III" = "#804FB3",
      "IV" = "#9969C7"
    ),
    na.value = NA
  )+
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()
ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_IUCN.png"), plot = p, dpi = 300)

p<-ggplot(grid%>%filter(!is.na(IUCN_CAT)),
          aes(x = IUCN_CAT,
              y = as.numeric(area))) +
  stat_summary(fun = sum, geom = "col", fill = "steelblue") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6)) +
  labs(
    x = "IUCN",
    y = expression(Area~(km^2))
  ) +
  theme_minimal() +
  theme(text = element_text(size = 20))

ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_area_stats_IUCN.png"), plot = p, width = 8, height = 6, dpi = 300)




# stats core prot per LULC
# lulc_core_prot<-grid%>%
#   group_by(sampled_habitat,core_prot_area_old)%>%
#   summarise(area_km2=sum(area)/10^6)%>%
#   st_drop_geometry()
# 
# lulc_gap<-grid%>%
#   group_by(sampled_habitat)%>%
#   summarise(tot_habitat_area = sum(as.numeric(area))/10^6,
#             target_core_prot_area = target_core_prot_fraction *tot_habitat_area)%>%
#   st_drop_geometry()
#   
# lulc_gap <- lulc_gap %>%
#   left_join(lulc_core_prot%>%filter(core_prot_area_old == T), by = "sampled_habitat")
# 
# lulc_gap$area_km2[is.na(lulc_gap$area_km2)] <- 0
# lulc_gap$gap_core_prot<-lulc_gap$target_core_prot_area - lulc_gap$area_km2
# 
# lulc_gap$rel_gap <- lulc_gap$gap_core_prot/lulc_gap$tot_habitat_area
write.csv(gap,paste0("WP4/2_output/01_PA_analysis/",siteID,"_gap_analysis.csv"))




### gap plot

gap_plot <- gap %>%
  dplyr::select(lulc_name, protected_area, gap_area) %>%
  tidyr::pivot_longer(
    cols = c(protected_area, gap_area),
    names_to = "status",
    values_to = "area"
  )

ggplot(gap_plot,
       aes(y = reorder(lulc_name, area),
           x = area/1e6,
           fill = status)) +
  
  geom_col(width = 0.8) +
  
  geom_errorbar(
    data = gap,
    aes(y = lulc_name,
        xmin = target_area/1e6,
        xmax = target_area/1e6),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  scale_fill_manual(values = c(
    protected_area = "#00A300",
    gap_area = "#D95F02"
  ),
  labels = c("Gap", "Existing core PA")) +
  
  labs(
    x = expression("Area (km"^2*")"),
    y = NULL,
    fill = NULL
  ) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )

ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_gap.png"), plot = p, width = 8, height = 6, dpi = 300)



# write grid for further processing
st_write(grid, paste0("WP4/2_output/02_optim/",siteID,"_input_grid.json"), driver = "GeoJSON", overwrite = T)



