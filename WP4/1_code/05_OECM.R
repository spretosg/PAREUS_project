### OECM suitability
library(ggpubr)
library(sf)
library(dplyr)
library(purrr)

# w_cult = 0.1
# w_prov =0.1
# w_ec = 0.5
# w_connect = 0.3

w_cult = 0.1
w_prov =0.1
w_ec = 0.4
w_connect = 0.4

source("WP4/1_code/wp4_functions_utils.R")
siteID<-"FRL04"
main_dir<-"P:/312204_pareus/"
stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/study_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="FRL04")


pu<-st_read(paste0("WP4/2_output/02_optim/PA_optim_conn_",siteID,".geojson"))

pu$oecm_suit<-oecm_lin_w(pu$sampled_cult_scaled,pu$sampled_prov_scaled,
                         pu$sampled_condition_scaled,pu$mw_connectivity,w_cult,w_prov,w_ec,w_connect)

# suit<-ggplot(pu) +
#   geom_sf(aes(fill = oecm_suit), color = NA) +
#   scale_fill_viridis_c(option = "viridis", name = "OECM suitability") +
#   geom_sf(data = stud_area, fill = NA, color = "black") +
#   theme_minimal()+
#   theme(text = element_text(size = 20))

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
             scen = "Base-lulc")
b<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_lulc",
                oecm_df=oecm_nat_lulc,
                stud_area,
                scen = "Nat-lulc")

c<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_global",
                oecm_df=oecm_nat_global,
                stud_area,
                scen = "Nat-glob")


d<-plot_pca_map(pu=pu,
                corePA_col = "core_pa_global",
                oecm_df=oecm_base_global,
                stud_area,
                scen = "base-glob")

library(ggnewscale)

# ggsave(paste0("WP4/2_output/03_pca_landscape/",siteID,"_optim_PCA_landscape.png"), plot = p, width = 12, height = 18, dpi = 300)
