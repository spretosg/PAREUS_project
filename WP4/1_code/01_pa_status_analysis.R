# GAP analysis
#how much of each natural land cover type is protected through IUCN cat Ia and II?
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)

# siteID<-"FRL04"
#siteID<-"SK021"
siteID<-"SK021"
target_core_prot_fraction <-0.1 #how much of each biome should be protected strictly
main_dir<-"P:/312204_pareus/"

#inputs

## stud_area 
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/stud_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="SK021")
target_crs <- 2154 #adjust this for the area
stud_area<-st_transform(stud_area,target_crs)
total_area<-st_area(stud_area)
stud_ara_vect <- vect(stud_area)

grid <- st_make_grid(stud_area, cellsize = 750, square = F)
grid <- st_sf(geometry = grid)


grid<-st_intersection(grid,stud_area["siteID"])
grid$ID<-c(1:nrow(grid))
grid$area<-as.numeric(st_area(grid))

# lulc
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

## PA
PA<-st_read(paste0(main_dir,"WP4/pa_existing/WDPA_SK021.shp"))
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


#sampling number of PA in one cell
# get intersection index list
ints <- st_intersects(grid, PA)

# count unique IUCN categories per grid cell
grid$n_pa <- sapply(ints, function(i) {
  length(unique(PA$IUCN_CAT[i]))
})

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

#sample highest status of protection
grid <- st_join(grid, PA[, c("class", "IUCN_CAT")], join = st_intersects)
grid <- grid %>%
  group_by(ID) %>%
  slice_max(class, n = 1, with_ties = FALSE) %>%
  ungroup()

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


### core prot vs other pa's
grid<-grid%>%mutate(core_prot_area_old = case_when(class<6 ~ F, class>5 ~ T, is.na(class)~F))

# stats core prot per LULC
lulc_core_prot<-grid%>%
  group_by(sampled_habitat,core_prot_area_old)%>%
  summarise(area_km2=sum(area)/10^6)%>%
  st_drop_geometry()

lulc_gap<-grid%>%
  group_by(sampled_habitat)%>%
  summarise(tot_habitat_area = sum(as.numeric(area))/10^6,
            target_core_prot_area = target_core_prot_fraction *tot_habitat_area)%>%
  st_drop_geometry()
  
lulc_gap <- lulc_gap %>%
  left_join(lulc_core_prot%>%filter(core_prot_area_old == T), by = "sampled_habitat")

lulc_gap$area_km2[is.na(lulc_gap$area_km2)] <- 0
lulc_gap$gap_core_prot<-lulc_gap$target_core_prot_area - lulc_gap$area_km2

lulc_gap$rel_gap <- lulc_gap$gap_core_prot/lulc_gap$tot_habitat_area
write.csv(lulc_gap,paste0("WP4/2_output/01_PA_analysis/",siteID,"_gap_analysis.csv"))




### gap plot

df_long <- lulc_gap %>%
  pivot_longer(
    cols = c(area_km2, gap_core_prot),
    names_to = "type",
    values_to = "area"
  )
df_long$type <- factor(df_long$type, levels = c("gap_core_prot", "area_km2"))

df_long<-df_long%>%mutate(lulc=case_when(sampled_habitat == 3 ~ "Forest",
                        sampled_habitat == 4 ~ "Wetlands",
                        sampled_habitat == 5 ~ "Water"))

p<-ggplot(df_long, aes(x = lulc, y = area, fill = type)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "gap_core_prot" = "#FDBE85",   # light orange,
      "area_km2" = "green"
      
    ),
    labels = c(
      "area_km2" = "current protected area",
      "gap_core_prot" = "Gap area"
    )
  ) +
  labs(
    x = "",
    y = "Area",
    fill = ""
  ) +
  theme_minimal()
ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_gap.png"), plot = p, width = 8, height = 6, dpi = 300)



# write grid for further processing
st_write(grid, paste0("WP4/2_output/02_optim/",siteID,"_input_grid.json"), driver = "GeoJSON", overwrite = T)



#### ---- Area statistics ----
#pa categories
# area_stats_IUCN<-PA%>%group_by(IUCN_CAT)%>%
#   summarise(area_km2 = sum(st_area(geometry)),
#             area_km2_uni = sum(st_area(st_union(geometry)))) %>%
#   mutate(area_km2 = as.numeric(area_km2) / 1e6, 
#          area_km2_uni = as.numeric(area_km2_uni) / 1e6)%>%st_drop_geometry()
# 
# 
# #plot

# 
# p<-ggplot() +
#   geom_sf(data = stud_area, fill = "grey95", color = "black") +
#   geom_sf(data = PA, aes(fill = IUCN_CAT), color = NA, alpha = 0.8) +
#   scale_fill_brewer(palette = "Set2", name = "IUCN Category") +
#   theme_minimal()
# ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_maps_IUCN.png"), plot = p, width = 8, height = 6, dpi = 300)
# 
# 
# 
# 
# 
# # increase lulc res to make better area estimates
# 
# 
# ## lulc classes area
# lulc_area <-  freq(lulc_highres)
# # area of one pixel (in map units)
# res_m <- res(lulc_highres)                 # pixel resolution
# cell_area <- res_m[1] * res_m[2]
# 
# # add area column
# lulc_area$lulc_tot_area_m2 <- lulc_area$count * cell_area
# lulc_area<-lulc_area%>%filter(value>0)
# 
# sum(lulc_area$lulc_tot_area_m2)-as.numeric(total_area)
# 
# 
# # LULC per IUCN poly
# results <- list()
# cell_area <- prod(res(lulc_highres))
# 
# for (cls in unique(iucn_union$strict_pa)) {
#   
#   poly_sub <- iucn_union[iucn_union$strict_pa == cls, ]
#   poly_sub <- st_collection_extract(poly_sub, "POLYGON")
#   sum(st_area(poly_sub))
#   
#   
#   r_mask <- mask(crop(lulc_highres, poly_sub), poly_sub)
#   
# 
#   area_tab <- expanse(r_mask, byValue = TRUE)
#   
#   area_tab$strict_pa <- cls
#   
#   results[[as.character(cls)]] <- area_tab
# }
# 
# final <- do.call(rbind, results)%>%filter(value>0)
# sum(final$area) - as.numeric(sum(st_area(iucn_union)))
# 
# # join total LULC val to final
# final<-merge(final,lulc_area,by="value")
# final$prot_frac<-final$area/final$lulc_tot_area_m2
# final<-final%>%mutate(LULC_class = case_when(value ==1 ~ "artificial_surfaces",
#                                              value ==2 ~ "agricultural_areas",
#                                              value ==3 ~ "forests_semi_natural",
#                                              value ==4 ~ "wetlands",
#                                              value ==5 ~ "water"))%>%dplyr::select(LULC_class,lulc_tot_area_m2,strict_pa,area,prot_frac)
# colnames(final)<-c("lulc_class","tot_lulc_area_m2","strict_pa","prot_area_m2","rel_protection")
# 
# gap_strict_pa<-final%>%filter(lulc_class %in% c("forests_semi_natural","water","wetlands") & strict_pa == T)%>%
#   mutate(
#     gap_strict_PA_km2 = ((tot_lulc_area_m2*0.1)-prot_area_m2)/10^6,
#     prot_area_km2 = prot_area_m2/10^6,
#     tot_lulc_area_km2 = tot_lulc_area_m2/10^6)%>%dplyr::select(lulc_class,prot_area_km2,gap_strict_PA_km2)
# 
# 
# 
# # p<-ggplot(gap_strict_pa,
# #        aes(x = reorder(lulc_class, -delta_strict_PA_km2),
# #            y = delta_strict_PA_km2)) +
# #   geom_col(fill = "steelblue") +
# #   labs(title = paste0(siteID," gap to protect of strict protected area (Ia/II) to protect 10% of the total LULC area"),
# #        x = "LULC",
# #        y = "km²") +
# #   theme_minimal(base_size = 14)
# df_long <- gap_strict_pa %>%
#   pivot_longer(cols = c(gap_strict_PA_km2, prot_area_km2),
#                names_to = "type",
#                values_to = "area")
# 
# p<-ggplot(df_long, aes(x = lulc_class, y = area, fill = type)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = c(
#     "gap_strict_PA_km2" = "red",
#     "prot_area_km2" = "green"
#   )) +
# #
#   theme_minimal()
# ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_gap_strict_prot.png"), plot = p, width = 8, height = 6, dpi = 300)
# 
# 
# 
# #### ---- Condition statistics 1####
# 
# 
# # 2. Extract ec values inside PA by IUCN category
# pa_vals <- extract(ec, PA)
# 
# # attach category
# pa_vals$IUCN_cat <- PA$IUCN_CAT[pa_vals$ID]
# 
# 
# # 3. Create raster of PA footprint
# pa_r <- rasterize(PA, ec)
# 
# # 4. Values outside PA
# ec_out <- mask(ec, pa_r, inverse=TRUE)
# 
# out_vals <- values(ec, na.rm=TRUE)
# out_df <- data.frame(
#   IUCN_cat = "Outside_PA",
#   ec = out_vals
# )
# colnames(out_df)[2] <- "ec"
# colnames(pa_vals)[2] <- "ec"
# 
# # 5. Combine
# df <- rbind(pa_vals%>%select(IUCN_cat,ec), out_df)
# 
# # 6. Boxplot
# p<-ggplot(df, aes(x=IUCN_cat, y=ec)) +
#   geom_boxplot() +
#   theme_bw() +
#   labs(x="IUCN Category", y="Ecosystem condition") +
#   theme(axis.text.x = element_text(angle=45, hjust=1))
# ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_EC_IUCN.png"), plot = p, width = 8, height = 6, dpi = 300)
# 
# #test
# set.seed(13)
# df_sample <- df[sample(nrow(df), 10000), ]
# kw<-kruskal.test(ec ~ IUCN_cat, data=df_sample)
# pair_test<-pairwise.wilcox.test(df_sample$ec, df_sample$IUCN_cat, p.adjust.method="BH")
# capture.output(
#   pair_test,
#   file = paste0("WP4/2_output/01_PA_analysis/",siteID,"_kw_pair.txt")
# )
