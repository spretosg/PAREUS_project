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

target_site<-"SK021"
####---- User parameter ----####
# weighting of OECM suitability raster components
# w_cult = 0.05
# w_prov =0.05
# w_ec = 0.45
# w_connect = 0.45
w_cult = 0.25
w_prov =0.25
w_ec = 0.25
w_connect = 0.25


coverage <- c(
  forests = 0.15,
  water = 0.2,
  wetland = 0.2,
  agricultural = 0.15
)

####---- Input and processing ----####
stud_area<-read_sf(paste0("data/shared/pareus_sites.gpkg"))%>%filter(siteID == target_site)

# PUs containing the optimized core PA
pu<-st_read(paste0("outputs/WP4/02_optim/PA_optim_grid_",target_site,".geojson"))

#connectivity raster (run 00_helper_connectivity.R first)
connectivity<-rast(paste0("outputs/WP4/02_optim/mw_connectivity_",target_site,".tif"))

# sample connetivity values to pu
pu$connectivity<- terra::extract(connectivity$normalized_current, pu, fun = mean, na.rm = TRUE)[,2]
pu<-zero_one_scale(
  pu,
  cols = c("connectivity")
)



####---- Linear combination for OECM suitability ----####
## check correlation of MCDA variables
#vars <- c("sampled_es_cult_scaled", "sampled_es_prov_scaled", "sampled_eco_cond_scaled", "connectivity_scaled")

# cor_mat <- pu %>%
#   st_drop_geometry() %>%
#   dplyr::select(dplyr::all_of(vars)) %>%
#   cor(use = "pairwise.complete.obs")
# 
# corrplot::corrplot(
#   cor_mat,
#   method = "color",
#   type = "upper",
#   addCoef.col = "black",   # show r values
#   number.cex = 0.8,        # size of r values
#   tl.col = "black",
#   tl.srt = 45,
#   diag = FALSE
# )

pu$oecm_suit<-oecm_lin_w(pu$sampled_es_cult_scaled,pu$sampled_es_prov_scaled,
                         pu$sampled_eco_cond_scaled,pu$connectivity_scaled,w_cult,w_prov,w_ec,w_connect)

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
oecm_nat_lulc<-select_oecm(pu=pu, 
                  mode= "class", 
                  coverage = coverage,
                  lulc_col = "lulc_name",
                  suitability_col = "oecm_suit",
                  corePA_col = "core_pa_lulc",
                  search_oecm_in = c("not protected","other protected areas"))

oecm_nat_global<-select_oecm(pu=pu, 
                           mode= "class", 
                           coverage = coverage,
                           lulc_col = "lulc_name",
                           suitability_col = "oecm_suit",
                           corePA_col = "core_pa_global",
                           search_oecm_in = c("not protected","other protected areas"))

oecm_base_lulc<-select_oecm(pu=pu, 
                  mode= "global", 
                  coverage = 0.2,
                  suitability_col = "oecm_suit",
                  corePA_col = "core_pa_lulc",
                  search_oecm_in = c("not protected","other protected areas"))

oecm_base_global<-select_oecm(pu=pu, 
                            mode= "global", 
                            coverage = 0.2,
                            suitability_col = "oecm_suit",
                            corePA_col = "core_pa_global",
                            search_oecm_in = c("not protected","other protected areas"))


####---- combine optimized PA and selected OECM into one PCA map, calc statistics ----####
base_lulc<-plot_pca_map(pu=pu,
             corePA_col = "core_pa_lulc",
             oecm_df=oecm_base_lulc,
             stud_area,
             scen = "BASE_LULC",
             save_output = T)
nat_lulc<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_lulc",
                oecm_df=oecm_nat_lulc,
                stud_area,
                scen = "NAT_LULC",
                save_output = T)

nat_glob<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_global",
                oecm_df=oecm_nat_global,
                stud_area,
                scen = "NAT_GLOB",
                save_output = T)


base_glob<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_global",
                oecm_df=oecm_base_global,
                stud_area,
                scen = "BASE_GLOB",
                save_output = T)
st_write(pu,paste0("outputs/WP4/03_pca_landscape/PCA_final",target_site,".geojson"))

stats_all<-rbind(base_lulc$stats,nat_lulc$stats,nat_glob$stats, base_glob$stats)
stats_all$pa_group[stats_all$pa_group == "Other PA (IUCN III-VI) high suitability for OECM"] <- "other_pa"
stats_all$pa_group[stats_all$pa_group == "other_pa"] <- "Other PA (IUCN III-VI)"
stats_all$pa_group[stats_all$pa_group == "core_PA"] <- "Optimized IUCN Ia or II"
stats_all$pa_group[stats_all$pa_group == "not_protected"] <- "No protection"
stats_all$pa_group[stats_all$pa_group == "High OECM suitability not protected - potential OECM"] <- "Potential OECM"
stats_all<-stats_all%>%group_by(lulc_name,scenario,pa_group)%>%summarise(area = sum(area)/10^6)


base_scenario <- "BASE_GLOB"

area_diff <- stats_all %>%
  group_by(lulc_name, pa_group) %>%
  mutate(
    base_area = area[scenario == base_scenario][1],
    rel_diff = 100 * (area - base_area) / base_area
  ) %>%
  ungroup()


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