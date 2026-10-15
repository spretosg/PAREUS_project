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
pu$oecm_suit<-oecm_lin_w(pu$sampled_es_cult_scaled,pu$sampled_es_prov_scaled,
                         pu$sampled_es_reg_scaled,pu$connectivity_scaled,human_c$w_cult,human_c$w_prov,human_c$w_reg,human_c$w_connect)


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


### SK LCP connectivity
#interpolate oecm suit if na
library(RANN)
pu <- fill_suit_na(pu, "oecm_suit", method = "neighbours")
sum(pu$suitability_filled_imputed) 

lcp <- lcp_connect(
  pu,
  sel            = pu$is_core_prot_agri,      # existing + new + upgrade cells
  suit_col       = "oecm_suit_filled",               # your suitability column
  exclude        = pu$lulc_name == "built-up",  # hard barrier
  transform      = "linear",
  max_resistance = 50,
  topology       = "mst"
)

pu$is_corridor_agri <- pu$id %in% lcp$corridor_ids

cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#bf8bff",
  "proposed upgrade existing PA" = "#e5d0ff",
  "not protected" = NA,
  "other protected areas" = "#F1EE8D"
)

ggplot() +
  geom_sf(data = pu[pu$is_corridor_agri, ], fill = "orange", colour = NA) +
  geom_sf(data = pu, aes(fill = core_pa_agri), color = NA) +
  scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  # geom_sf(data = lcp$paths, colour = "red", linewidth = 0.6) +
  geom_sf(data = lcp$nodes, size = 1.5) +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  geom_sf(
    data = pu%>%filter(lulc_name == "built-up"),
    aes(fill = lulc_name),
    color = "grey", fill="grey"  )+
  #scale_fill_manual(values = cols, breaks = names(cols)[1:3], name = NULL) +
  theme_minimal()



#path statistics
corridor_lulc_stats <- function(pu, flag_col, scenario) {
  pu %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(.data[[flag_col]] %in% TRUE) %>%
    dplyr::group_by(lulc_name) %>%
    dplyr::summarise(SQkm = sum(area_km2), n_cells = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(pct = 100 * SQkm / sum(SQkm), scenario = scenario)
}

# adjust the labels if your naming is the other way round
sc_levels <- c("No PA's in agricultural areas", "Agricultural areas as PA possible")

stats <- dplyr::bind_rows(
  corridor_lulc_stats(pu, "is_corridor_agri", sc_levels[1]),
  corridor_lulc_stats(pu, "is_corridor_global",   sc_levels[2])
)

comparison <- stats %>%
  dplyr::select(lulc_name, scenario, SQkm) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = SQkm, values_fill = 0)
comparison <- dplyr::bind_rows(
  comparison,
  comparison %>% dplyr::summarise(dplyr::across(where(is.numeric), sum)) %>%
    dplyr::mutate(lulc_name = "TOTAL")
)
comparison

lulc_order <- stats %>%
  dplyr::group_by(lulc_name) %>% dplyr::summarise(t = sum(SQkm)) %>%
  dplyr::arrange(dplyr::desc(t)) %>% dplyr::pull(lulc_name)

pies <- stats %>%
  dplyr::mutate(lulc_name = factor(lulc_name, levels = lulc_order),
                scenario  = factor(scenario,  levels = sc_levels)) %>%
  dplyr::arrange(scenario, lulc_name) %>%
  dplyr::group_by(scenario) %>%
  dplyr::mutate(total = sum(SQkm),
                frac  = SQkm / total,
                end   = 2 * pi * cumsum(frac),
                start = end - 2 * pi * frac,
                mid   = (start + end) / 2) %>%
  dplyr::ungroup()

# pie area proportional to total corridor area; pies placed side by side
gap <- 0.8
sizes <- pies %>%
  dplyr::distinct(scenario, total) %>%
  dplyr::arrange(scenario) %>%
  dplyr::mutate(r  = sqrt(total / max(total)),
                x0 = cumsum(2 * r + gap) - r - gap)

pies <- dplyr::left_join(pies, sizes[, c("scenario", "r", "x0")], by = "scenario")

min_label_pct <- 1   # skip labels for slices smaller than this
labs_df <- pies %>%
  dplyr::filter(frac >= min_label_pct / 100) %>%
  dplyr::mutate(lx  = x0 + 1.12 * r * sin(mid),
                ly  = 1.12 * r * cos(mid),
                lab = paste0(format(round(SQkm, 1), nsmall = 1), " km²"),
                hj  = ifelse(sin(mid) >= 0, 0, 1))

titles <- sizes %>%
  dplyr::mutate(lab = sprintf("%s\n%.0f km² total", scenario, total),
                ty  = -r - 0.3)

# 3. Plot -------------------------------------------------------------------
map_pie <- ggplot() +
  ggforce::geom_arc_bar(
    data = pies,
    aes(x0 = x0, y0 = 0, r0 = 0, r = r, start = start, end = end,
        fill = lulc_name, group = interaction(scenario, lulc_name)),
    colour = "white", linewidth = 0.3
  ) +
  geom_text(data = labs_df, aes(lx, ly, label = lab, hjust = hj), size = 3.5) +
  geom_text(data = titles, aes(x0, ty, label = lab), fontface = "bold", lineheight = 0.9) +
  scale_fill_viridis_d(name = "LULC", drop = FALSE) +   # swap for your own LULC palette
  coord_fixed(clip = "off") +
  labs(subtitle = "Pie area proportional to total corridor area") +
  theme_void(base_size = 14) +
  theme(legend.position = "bottom",
        plot.margin = margin(10, 40, 10, 40))

map_pie

## corridor growth where oecm suitability is best to meet a target


### FR and NOR
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