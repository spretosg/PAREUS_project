############################################################################
# Pareus WP4 - OECM and PCA landscape
# Date: 06.08.2026
# Author: Reto Spielhofer (Norwegian Institute for Nature Research)
############################################################################
library(ggpubr)
library(sf)
library(dplyr)
library(purrr)
library(terra)
library(ggnewscale)
library(tidyverse)

source("code/WP4/wp4_functions_utils.R")

target_site<-"TRD"
scenario<-"nature"# "nature" #"intermediate" "multifunctional"
####---- User parameter ----####
#MCDA weighting 
human_c<-list(
  w_cult = 0.4,
  w_prov =0.4,
  w_reg = 0.1,
  w_connect = 0.1
)

eco_c<-list(
  w_cult = 0.1,
  w_prov =0.1,
  w_reg = 0.4,
  w_connect = 0.4
)


coverage <- c(
  forest = 0.2,
  water = 0.2,
  wetland = 0.2,
  agricultural = 0.2
)

####---- Input and processing ----####
stud_area<-read_sf(paste0("data/shared/",target_site,".gpkg"))

# PUs containing the optimized core PA
pu<-st_read(paste0("outputs/WP4/02_optim/PA_optim_grid_",target_site,"_",scenario,".geojson"))
lulc_stats<-pu%>%st_drop_geometry()%>%
  group_by(lulc_name)%>%summarize(sum_area = sum(area_km2))

#connectivity raster (run 00_helper_connectivity.R first)
connectivity<-rast(paste0("outputs/WP4/02_optim/conn_",target_site,".tif"))

# sample connetivity values to pu
pu$connectivity<- terra::extract(connectivity$normalized_current, pu, fun = mean, na.rm = TRUE)[,2]
pu<-zero_one_scale(
  pu,
  cols = c("connectivity")
)

#eco centric
pu$oecm_suit<-oecm_lin_w(pu$sampled_es_cult_scaled,pu$sampled_es_prov_scaled,
                         pu$sampled_es_reg_scaled,pu$connectivity_scaled,eco_c$w_cult,eco_c$w_prov,eco_c$w_reg,eco_c$w_connect)

#human centric
# pu$oecm_suit<-oecm_lin_w(pu$sampled_es_cult_scaled,pu$sampled_es_prov_scaled,
#                          pu$sampled_es_reg_scaled,pu$connectivity_scaled,human_c$w_cult,human_c$w_prov,human_c$w_reg,human_c$w_connect)


ggplot(pu %>% filter(!existing_corePA)) +
  geom_sf(aes(fill = oecm_suit), color = NA) +
  scale_fill_gradientn(
    colours = c(
      "#ADD8E6",  # light blue
      "#388ca7",  # blue
      "#6A3D9A",  # purple
      "#4B0082"   # dark purple
    ),
    limits = c(0, 1),
    name = "OECM suitability")+
  #scale_fill_viridis_c(option = "viridis", name = "OECM suitability") +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 15),
        legend.position = "top")


####---- Select top-% OECM suitability based on scenario----####
# oecm_nat_lulc<-select_oecm(pu=pu, 
#                   mode= "class", 
#                   coverage = coverage,
#                   lulc_col = "lulc_name",
#                   suitability_col = "oecm_suit",
#                   corePA_col = "core_pa_lulc",
#                   search_oecm_in = c("not protected","other protected areas"))

# oecm_nat_global<-select_oecm(pu=pu, 
#                            mode= "class", 
#                            coverage = coverage,
#                            lulc_col = "lulc_name",
#                            suitability_col = "oecm_suit",
#                            corePA_col = "core_pa_global",
#                            search_oecm_in = c("not protected","other protected areas"))

oecm_lulc<-select_oecm(pu=pu, 
                  mode= "class", 
                  coverage = coverage,
                  lulc_col = "lulc_name",
                  suitability_col = "oecm_suit",
                  corePA_col = "core_pa_lulc",
                  search_oecm_in = c("not protected","other protected areas"),
                  w= 0.9, order = 2)

oecm_global<-select_oecm(pu=pu, 
                            mode= "global", 
                            coverage = 0.2,
                            suitability_col = "oecm_suit",
                            corePA_col = "core_pa_global",
                            search_oecm_in = c("not protected","other protected areas"),
                         w= 0.9, order = 2)


####---- combine optimized PA and selected OECM into one PCA map, calc statistics ----####
base_lulc<-plot_pca_map(pu=pu,
                        corePA_col = "core_pa_lulc",
                        oecm_df=oecm_lulc,
                        stud_area,
                        scen = "LULC PA/OECM",
                        save_output = F)

base_glob<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_global",
                oecm_df=oecm_global,
                stud_area,
                scen = "GLOB PA/OECM",
                save_output = F)
st_write(pu,paste0("outputs/WP4/03_pca_landscape/PCA_final",target_site,".geojson"))

stats_all<-rbind(base_lulc$stats,base_glob$stats)
stats_all$pca[stats_all$pca == "Other PA (IUCN III-VI) high suitability for OECM"] <- "other_pa"
stats_all$pca[stats_all$pca == "other_pa"] <- "Other PA (IUCN III-VI)"
stats_all$pca[stats_all$pca == "core_PA"] <- "Optimized IUCN Ia or II"
stats_all$pca[stats_all$pca == "not_protected"] <- "No protection"
stats_all$pca[stats_all$pca == "High OECM suitability not protected - potential OECM"] <- "Potential OECM"
stats_all<-stats_all%>%group_by(lulc_name,scenario,pca)%>%summarise(area = sum(area)/10^6)


base_scenario <- "GLOB PA/OECM"

area_diff <- stats_all %>%
  group_by(lulc_name, pca) %>%
  mutate(
    base_area = area[scenario == base_scenario][1],
    rel_diff = 100 * (area - base_area) / base_area
  ) %>%
  ungroup()%>%filter(scenario != base_scenario)




# ggplot(area_diff %>% filter(!lulc_name %in% c(NA,"built-up") & scenario != "BASE_GLOB"),
#        aes(x = scenario,
#            y = pa_group,
#            fill = rel_diff)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = sprintf("%.0f%%", rel_diff)), size = 3) +
#   facet_wrap(~lulc_name) +
#   scale_fill_gradient2(
#     low = "#B2182B",
#     mid = "white",
#     high = "#2166AC",
#     midpoint = 0,
#     name = "Change (%)"
#   ) +
#   theme_bw()