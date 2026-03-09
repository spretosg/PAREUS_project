# GAP analysis
#how much of each natural land cover type is protected through IUCN cat Ia and II?
library(sf)
library(terra)
library(dplyr)
library(ggplot2)

# siteID<-"FRL04"
siteID<-"FRA_BAR2"

main_dir<-"P:/312204_pareus/"

#inputs

## stud_area 
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/stud_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="FRA_BAR2")
target_crs <- 2154 #adjust this for the area
stud_area<-st_transform(stud_area,target_crs)
total_area<-st_area(stud_area)
stud_ara_vect <- vect(stud_area)

# lulc
lulc<-paste0(main_dir,"WP4/habitat/",siteID,"_lulc.tif")
lulc<-terra::rast(lulc)
lulc <- project(lulc, paste0("epsg:",target_crs))
lulc <- crop(lulc, stud_ara_vect)
lulc <- mask(lulc, stud_ara_vect)

# EC for differences in EC between PA
ec<-paste0(main_dir,"WP4/features/",siteID,"_ec.tif")
ec<-terra::rast(ec)
ec <- project(ec, paste0("epsg:",target_crs))
ec <- crop(ec, stud_ara_vect)
ec <- mask(ec, stud_ara_vect)

# PA
PA<-st_read(paste0(main_dir,"WP4/pa_existing/WDPA_",siteID,".shp"))
PA<-st_transform(PA,target_crs)
PA<-st_make_valid(PA)

#crop to study area
PA<-st_intersection(PA,stud_area)

#### ---- Area statistics ----
#pa categories
area_stats_IUCN<-PA%>%group_by(IUCN_CAT)%>%
  summarise(area_km2 = sum(st_area(geometry)),
            area_km2_uni = sum(st_area(st_union(geometry)))) %>%
  mutate(area_km2 = as.numeric(area_km2) / 1e6, 
         area_km2_uni = as.numeric(area_km2_uni) / 1e6)%>%st_drop_geometry()


#plot
p<-ggplot(area_stats_IUCN, 
       aes(x = reorder(IUCN_CAT, -area_km2_uni),
           y = area_km2_uni)) +
  geom_col(fill = "steelblue") +
  labs(title = paste0(siteID," Area by IUCN Category"),
        x = "IUCN Category",
       y = "Area (km²)") +
  theme_minimal(base_size = 14)
  
ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_area_stats_IUCN.png"), plot = p, width = 8, height = 6, dpi = 300)

p<-ggplot() +
  geom_sf(data = stud_area, fill = "grey95", color = "black") +
  geom_sf(data = PA, aes(fill = IUCN_CAT), color = NA, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", name = "IUCN Category") +
  theme_minimal()
ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_maps_IUCN.png"), plot = p, width = 8, height = 6, dpi = 300)



###reclass lulc
lulc_recl<- floor(lulc / 100)
#plot(lulc_recl)

# reclass PA importance
PA <- PA %>%
  mutate(num_degree_prot = case_when(
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



split_PA <- st_intersection(PA) %>%
  group_by(geometry) %>%
  slice_max(num_degree_prot, n = 1) %>%
  ungroup()%>%dplyr::select(IUCN_CAT,num_degree_prot)

sum(st_area(split_PA))
iucn_union <- split_PA %>% mutate(strict_pa = case_when(num_degree_prot <6 ~ F,
                                                        num_degree_prot>5 ~ T))%>%
  group_by(strict_pa) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")
sum(st_area(iucn_union))

# increase lulc res to make better area estimates
lulc_highres <- disagg(lulc_recl, fact = 5, method = "near")

## lulc classes area
lulc_area <-  freq(lulc_highres)
# area of one pixel (in map units)
res_m <- res(lulc_highres)                 # pixel resolution
cell_area <- res_m[1] * res_m[2]

# add area column
lulc_area$lulc_tot_area_m2 <- lulc_area$count * cell_area
lulc_area<-lulc_area%>%filter(value>0)

sum(lulc_area$lulc_tot_area_m2)-as.numeric(total_area)


# LULC per IUCN poly
results <- list()
cell_area <- prod(res(lulc_highres))

for (cls in unique(iucn_union$strict_pa)) {
  
  poly_sub <- iucn_union[iucn_union$strict_pa == cls, ]
  poly_sub <- st_collection_extract(poly_sub, "POLYGON")
  sum(st_area(poly_sub))
  
  
  r_mask <- mask(crop(lulc_highres, poly_sub), poly_sub)
  

  area_tab <- expanse(r_mask, byValue = TRUE)
  
  area_tab$strict_pa <- cls
  
  results[[as.character(cls)]] <- area_tab
}

final <- do.call(rbind, results)%>%filter(value>0)
sum(final$area) #~app 40km2 off to the sum ~1.3%

# join total LULC val to final
final<-merge(final,lulc_area,by="value")
final$prot_frac<-final$area/final$lulc_tot_area_m2
final<-final%>%mutate(LULC_class = case_when(value ==1 ~ "artificial_surfaces",
                                             value ==2 ~ "agricultural_areas",
                                             value ==3 ~ "forests_semi_natural",
                                             value ==4 ~ "wetlands",
                                             value ==5 ~ "water"))%>%dplyr::select(LULC_class,lulc_tot_area_m2,strict_pa,area,prot_frac)
colnames(final)<-c("lulc_class","tot_lulc_area_m2","strict_pa","prot_area_m2","rel_protection")

gap_strict_pa<-final%>%filter(lulc_class %in% c("forests_semi_natural","water","wetlands") & strict_pa == T)%>%
  mutate(delta_strict_PA_km2 = ((tot_lulc_area_m2*0.1)-(tot_lulc_area_m2*rel_protection))/10^6)

p<-ggplot(gap_strict_pa, 
       aes(x = reorder(lulc_class, -delta_strict_PA_km2),
           y = delta_strict_PA_km2)) +
  geom_col(fill = "steelblue") +
  labs(title = paste0(siteID," gap to protect of strict protected area (Ia/II) to protect 10% of the total LULC area"),
       x = "LULC",
       y = "km²") +
  theme_minimal(base_size = 14)
ggsave(paste0("WP4/2_output/01_PA_analysis/",siteID,"_gap_strict_prot.png"), plot = p, width = 8, height = 6, dpi = 300)


write.csv(final,paste0(siteID,"_gap_analysis.csv"))

#### ---- Condition statistics ####


# 2. Extract ec values inside PA by IUCN category
pa_vals <- extract(ec, PA)

# attach category
pa_vals$IUCN_cat <- PA$IUCN_CAT[pa_vals$ID]

# remove ID column
pa_vals <- pa_vals[, c("IUCN_CAT", names(int))]

# 3. Create raster of PA footprint
pa_r <- rasterize(PA, ec)

# 4. Values outside PA
ec_out <- mask(ec, pa_r, inverse=TRUE)

out_vals <- values(ec, na.rm=TRUE)
out_df <- data.frame(
  IUCN_cat = "Outside_PA",
  ec = out_vals
)
colnames(out_df)[2] <- "ec"

# 5. Combine
colnames(pa_vals)[2] <- "ec"
df <- rbind(pa_vals%>%select(IUCN_cat,ec), out_df)

# 6. Boxplot
ggplot(df, aes(x=IUCN_cat, y=ec)) +
  geom_boxplot() +
  theme_bw() +
  labs(x="IUCN Category", y="EC value") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

#test
set.seed(1)
df_sample <- df[sample(nrow(df), 5000), ]
kruskal.test(ec ~ IUCN_cat, data=df_sample)
pairwise.wilcox.test(df_sample$ec, df_sample$IUCN_cat, p.adjust.method="BH")
