### OECM suitability
library(ggpubr)
library(sf)
library(dplyr)
library(purrr)
library(ggnewscale)


w_cult = 0.05
w_prov =0.05
w_ec = 0.45
w_connect = 0.45

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
pu<-zero_one_scale(
  pu,
  cols = c("mw_connectivity")
)
## check correlation of MCDA variables
vars <- c("sampled_cult_scaled", "sampled_prov_scaled", "sampled_condition_scaled", "mw_connectivity_scaled")

cor_mat <- pu %>%
  st_drop_geometry() %>%
  dplyr::select(dplyr::all_of(vars)) %>%
  cor(use = "pairwise.complete.obs")

corrplot::corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  addCoef.col = "black",   # show r values
  number.cex = 0.8,        # size of r values
  tl.col = "black",
  tl.srt = 45,
  diag = FALSE
)

pu$oecm_suit<-oecm_lin_w(pu$sampled_cult_scaled,pu$sampled_prov_scaled,
                         pu$sampled_condition_scaled,pu$mw_connectivity_scaled,w_cult,w_prov,w_ec,w_connect)

ggplot(pu %>% filter(!core_prot_area_old)) +
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

ggplot(pu, aes(x = lulc_name, y = oecm_suit, fill = lulc_name, group = lulc_name)) +
  geom_boxplot() +
  # scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  theme_minimal() +
  theme(legend.position = "none",text = element_text(size = 20))+
  labs(
    x = "",
    y = "OECM suitability",
    fill = ""
  )

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
# conn<-ggplot(pu) +
#   geom_sf(aes(fill = mw_connectivity), color = NA) +
#   scale_fill_viridis_c(option = "plasma", name = "Connectivity") +
#   geom_sf(data = stud_area, fill = NA, color = "black") +
#   theme_minimal()+
#   theme(text = element_text(size = 20))

# maps<-ggarrange(conn, suit)
# ggsave(paste0("WP4/2_output/02_optim/",siteID,"_OECM_connect.png"), plot = maps, width = 18, height = 10, dpi = 300)
# 

#### other plots
# ec<-ggplot(pu) +
#   geom_sf(aes(fill = sampled_condition_scaled), color = NA) +
#   scale_fill_viridis_c(option = "greens", name = "ecosystem condition") +
#   geom_sf(data = stud_area, fill = NA, color = "black") +
#   theme_minimal()+
#   theme(text = element_text(size = 20))
# ec
# 
# 
# cost_pol<-ggplot(pu) +
#   geom_sf(aes(fill = sampled_cost_pol), color = NA) +
#   scale_fill_viridis_c(option = "greens", name = "policy costs") +
#   geom_sf(data = stud_area, fill = NA, color = "black") +
#   theme_minimal()+
#   theme(text = element_text(size = 20))
# cost_pol
# 
# 
# ec<-ggplot(pu) +
#   geom_sf(aes(fill = sampled_condition_scaled), color = NA) +
#   scale_fill_viridis_c(option = "greens", name = "ecosystem condition") +
#   geom_sf(data = stud_area, fill = NA, color = "black") +
#   theme_minimal()+
#   theme(text = element_text(size = 20))
# ec



#exclude existing and new core pa here
# oecm <- pu %>% filter(is_core_prot_lulc == F & oecm_suit >= quantile(oecm_suit, 1-p_oecm, na.rm = TRUE))

coverage <- c(
  forests = 0.15,
  water = 0.2,
  wetland = 0.2,
  agricultural = 0.15
)

oecm_nat_lulc<-select_oecm(pu=pu, 
                  mode= "class", 
                  coverage = coverage,
                  lulc_col = "lulc_name",
                  suitability_col = "oecm_suit",
                  corePA_col = "core_pa_lulc",
                  exclude_corePA = T)

oecm_nat_global<-select_oecm(pu=pu, 
                           mode= "class", 
                           coverage = coverage,
                           lulc_col = "lulc_name",
                           suitability_col = "oecm_suit",
                           corePA_col = "core_pa_global",
                           exclude_corePA = T)

oecm_base_lulc<-select_oecm(pu=pu, 
                  mode= "global", 
                  coverage = 0.2,
                  suitability_col = "oecm_suit",
                  corePA_col = "core_pa_lulc",
                  exclude_corePA = T)

oecm_base_global<-select_oecm(pu=pu, 
                            mode= "global", 
                            coverage = 0.2,
                            suitability_col = "oecm_suit",
                            corePA_col = "core_pa_global",
                            exclude_corePA = T)


a<-plot_pca_map(pu=pu,
             corePA_col = "core_pa_lulc",
             oecm_df=oecm_base_lulc,
             stud_area,
             scen = "BASE_LULC")
b<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_lulc",
                oecm_df=oecm_nat_lulc,
                stud_area,
                scen = "NAT_LULC")

c<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_global",
                oecm_df=oecm_nat_global,
                stud_area,
                scen = "NAT_GLOB")


d<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_global",
                oecm_df=oecm_base_global,
                stud_area,
                scen = "BASE_GLOB")

stats_all<-rbind(a$stats,b$stats,c$stats, d$stats)
stats_all$pa_group[stats_all$pa_group == "Other PA (IUCN III-VI) high suitability for OECM"] <- "other_pa"
stats_all$pa_group[stats_all$pa_group == "other_pa"] <- "Other PA (IUCN III-VI)"
stats_all$pa_group[stats_all$pa_group == "core_PA"] <- "Optimized IUCN Ia or II"
stats_all$pa_group[stats_all$pa_group == "not_protected"] <- "No protection"
stats_all$pa_group[stats_all$pa_group == "High OECM suitability not protected - potential OECM"] <- "Potential OECM"
stats_all<-stats_all%>%group_by(lulc_name,scenario,pa_group)%>%summarise(area = sum(area))

library(dplyr)
library(tidyr)

base_scenario <- "BASE_GLOB"

area_diff <- stats_all %>%
  group_by(lulc_name, pa_group) %>%
  mutate(
    base_area = area[scenario == base_scenario][1],
    rel_diff = 100 * (area - base_area) / base_area
  ) %>%
  ungroup()


ggplot(area_diff %>% filter(!lulc_name %in% c(NA,"built-up") & scenario != "BASE_GLOB"),
       aes(x = scenario,
           y = pa_group,
           fill = rel_diff)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", rel_diff)), size = 3) +
  facet_wrap(~lulc_name) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    name = "Change (%)"
  ) +
  theme_bw()

library(ggnewscale)

# ggsave(paste0("WP4/2_output/03_pca_landscape/",siteID,"_optim_PCA_landscape.png"), plot = p, width = 12, height = 18, dpi = 300)
