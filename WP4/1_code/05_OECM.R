### OECM suitability
library(ggpubr)
w_cult = 0.25
w_prov =0.25
w_ec = 0.25
w_connect = 0.25
source("WP4/1_code/wp4_functions_utils.R")
siteID<-"FRL04"
main_dir<-"P:/312204_pareus/"
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/study_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="FRL04")


pu<-st_read(paste0("WP4/2_output/02_optim/PA_optim_conn_",siteID,".geojson"))

pu$oecm_suit<-oecm_lin_w(pu$sampled_cult_scaled,pu$sampled_prov_scaled,
                         pu$sampled_condition_scaled,pu$mw_connectivity,w_cult,w_prov,w_ec,w_connect)

suit<-ggplot(pu) +
  geom_sf(aes(fill = oecm_suit), color = NA) +
  scale_fill_viridis_c(option = "viridis", name = "OECM suitability") +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 20))

# ggplot(pu, aes(x = LULC_class, y = oecm_suit, fill = LULC_class, group = LULC_class)) +
#   geom_boxplot() +
#   # scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
#   theme_minimal() +
#   theme(legend.position = "none",text = element_text(size = 20))+
#   labs(
#     x = "",
#     y = "OECM suitability",
#     fill = ""
#   )

# conn<-ggplot(pu, aes(x = LULC_class, y = mw_connectivity, fill = LULC_class, group = LULC_class)) +
#   geom_boxplot() +
#   # scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
#   theme_minimal() +
#   theme(legend.position = "none",text = element_text(size = 20))+
#   labs(
#     x = "",
#     y = "Connectivity",
#     fill = ""
#   )
conn<-ggplot(pu) +
  geom_sf(aes(fill = mw_connectivity), color = NA) +
  scale_fill_viridis_c(option = "plasma", name = "Connectivity") +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 20))

maps<-ggarrange(conn, suit)
ggsave(paste0("WP4/2_output/02_optim/",siteID,"_OECM_connect.png"), plot = maps, width = 18, height = 10, dpi = 300)


#### other plots
ec<-ggplot(pu) +
  geom_sf(aes(fill = sampled_condition_scaled), color = NA) +
  scale_fill_viridis_c(option = "greens", name = "ecosystem condition") +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 20))
ec


cost_pol<-ggplot(pu) +
  geom_sf(aes(fill = sampled_cost_pol), color = NA) +
  scale_fill_viridis_c(option = "greens", name = "policy costs") +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 20))
cost_pol


ec<-ggplot(pu) +
  geom_sf(aes(fill = sampled_condition_scaled), color = NA) +
  scale_fill_viridis_c(option = "greens", name = "ecosystem condition") +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 20))
ec

p_oecm<-0.2 # from total area

#exclude existing and new core pa here
oecm <- pu %>% filter(is_core_prot_lulc == F & oecm_suit >= quantile(oecm_suit, 1-p_oecm, na.rm = TRUE))

# oecm <- oecm %>% filter(is_core_prot_lulc == F)
core_pa <- pu %>% filter(is_core_prot_lulc == T)

other_pa<-pu%>%filter(class<6)


## highest 20% OECM suitability

oecm$pa_group <- ifelse(oecm$n_pa > 0, "High OECM suitability already protected (IUCN III-VI)", "High OECM suitability not protected - potential OECM")

## other pa, remove ids from OECM
other_pa_filt <- other_pa[!other_pa$ID %in% oecm$ID, ]

other_pa_filt <- other_pa_filt[!other_pa_filt$ID %in% core_pa$ID, ]
other<-pu %>% filter(core_pa_lulc == "other")
other <- other[!other$ID %in% oecm$ID, ]
other <- other[!other$ID %in% other_pa_filt$ID, ]

stats_oecm<-oecm%>%group_by(pa_group,sampled_habitat)%>%summarise(area = sum(area))%>%st_drop_geometry()
stats_other_pa<-other_pa_filt%>%group_by(sampled_habitat)%>%summarise(area = sum(area))%>%st_drop_geometry()
stats_core_pa<-core_pa%>%group_by(sampled_habitat)%>%summarise(area = sum(area))%>%st_drop_geometry()
stats_other<-other%>%group_by(sampled_habitat)%>%summarise(area = sum(area))%>%st_drop_geometry()

total<-pu%>%group_by(sampled_habitat)%>%summarise(area = sum(area))%>%st_drop_geometry()

library(ggnewscale)
cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#00FF00",
  "proposed upgrade existing PA" = "#B8FFB8"
)

# plot(st_geometry(oecm),col="purple",add = T)
# plot(st_geometry(core_pa),col="red",add=T)
# 
# plot(st_geometry(other_pa_filt),col="lightblue",add=T)
# plot(st_geometry(other),col="lightgrey")

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
