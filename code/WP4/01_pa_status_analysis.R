# GAP analysis
#how much of each natural land cover type is protected through IUCN cat Ia and II?
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)
source("WP4/1_code/wp4_functions_utils.R")
siteID<-"FRL04"
main_dir<-"P:/312204_pareus/"

####---- User parameter ----####
#Parameters to set, ev to develop in collaboration with stakeholders
core_PA<-c("Ia","II") #WDPA categories to define as strictly protected areas

#the targets depends on the scenario, a global scenario defines a global protection target across all LULC, the lulc targets are used if lulc specific protection is aimed
target_glob <- 0.1 #global target for core protection across all LULC types
target_lulc <-c(forest = 0.1, water = 0.2, wetland = 0.2, agricultural = 0.05) #lulc specific protection target
PU_size<-1200 #lin m2 lower number increase the resolution but increase also computation time down stream

####---- Input and processing ----####
# study area
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/study_site.gpkg"))
# LULC raster
lulc<-terra::rast(paste0(main_dir,"WP4/habitat/",siteID,"_lulc.tif"))
# current PA network
PA<-st_read(paste0(main_dir,"WP4/pa_existing/WDPA_",siteID,".shp"))

## spatial transformation
target_crs<-25833
stud_area<-st_transform(stud_area,target_crs)
lulc <- project(lulc, paste0("epsg:",target_crs))
PA<-st_transform(PA,st_crs(target_crs))

#crop to study area
stud_area<-stud_area%>%filter(siteID=="FRL04")
lulc <- crop(lulc, vect(stud_area))
lulc <- mask(lulc, vect(stud_area))
PA <- PA%>%st_intersection(stud_area)%>%st_make_valid()

## establish a grid of planning units PU with an ID
grid <- stud_area %>%
  st_make_grid(cellsize = PU_size, square = FALSE) %>%
  st_sf(geometry = ., id = seq_along(.)) %>%
  st_intersection(stud_area["siteID"])

grid$area<-as.numeric(st_area(grid))

# extract habitat per PU 
#only use LULC main classes
lulc <- lulc %>%
  (\(x) floor(x / 100))() %>%
  disagg(fact = 5, method = "near")

grid <- grid %>%
  mutate(
    sampled_habitat = terra::extract(
      lulc,
      terra::vect(st_centroid(.))
    )[, 2]
  )%>%filter(!is.na(sampled_habitat) & sampled_habitat != 0)

## LULC statistics
lulc_stats<-grid%>%st_drop_geometry()%>%group_by(sampled_habitat)%>%summarise(area_km2=sum(as.numeric(area))/10^6)

## sample PA status 
PA <- PA %>%
  mutate(prot_class = case_when(
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


grid<-grid%>%mutate(lulc_name = case_when(sampled_habitat == 1 ~ "built-up",
                                          sampled_habitat == 2 ~ "agricultural",
                                          sampled_habitat == 3 ~ "forest",
                                          sampled_habitat == 4 ~ "wetland",
                                          sampled_habitat == 5 ~"water"))

#sampling number of PA in one cell
# get intersection index list
ints <- st_intersects(grid, PA)

grid <- grid %>%
  mutate(
    n_pa = sapply(ints, \(i) length(unique(PA$IUCN_CAT[i])))
  ) %>%
  st_join(PA[, c("prot_class", "IUCN_CAT")], join = st_intersects) %>%
  group_by(id) %>%
  slice_max(prot_class, n = 1, with_ties = FALSE) %>%
  ungroup()

### core prot vs other pa's
grid <- grid %>%
  mutate(
    existing_corePA = IUCN_CAT %in% core_PA
  )
# write grid for further processing
st_write(grid, paste0("WP4/2_output/01_PA_analysis/",siteID,"_input_grid.json"), driver = "GeoJSON", overwrite = T)

####---- Protection analysis ----####
p<-ggplot() +
     geom_sf(
       data = grid%>%filter(n_pa>0),
        aes(fill = n_pa),
     color = "NA"  )+
  scale_fill_gradient(
    high = "#C50202",
    low = "#FFC2C2",
    space = "Lab",
    na.value = "grey50",
  )+
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 20))
ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_n_pa.png"), plot = p, width = 8, height = 6, dpi = 300)

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

## area protected class
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


####---- Protection gap analysis ----####

gap_glob <- protection_gap(
  grid,
  lulc_col = "lulc_name",
  corePA_col = "exisiting_corePA",
  lockout = "built-up",
  mode = "global",
  target = target_glob
)

gap_lulc <- protection_gap(
  grid,
  lulc_col = "lulc_name",
  corePA_col = "exisiting_corePA",
  lockout = "built-up",
  mode = "class",
  target = target_lulc
)

gap<-rbind(gap_lulc,gap_glob)%>%filter(!is.na(lulc_name))
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