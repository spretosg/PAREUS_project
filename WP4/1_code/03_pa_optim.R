# test prioritizr
library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(gdistance)

main_dir<-"P:/312204_pareus/"
siteID<-"FRL04"
gap<-read.csv(paste0("WP4/2_output/01_PA_analysis/",siteID,"_gap_analysis.csv"))

stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/study_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="FRL04")

pu<-st_read(paste0("WP4/2_output/02_optim/",siteID,"_input_final_grid.json"))
pu$ID<-c(1:nrow(pu))
pu$area<-as.numeric(st_area(pu))
pu$inv_cost_pol<-1/pu$sampled_cost_pol

# pu<-pu%>%mutate(lock_out = case_when(
#   is.na(class) ~  FALSE,
#   class>5 ~ TRUE,
#   class<6 ~ FALSE
# ))


pu<-pu%>%mutate(inv_dist = case_when(
  min_distance_scaled == 0 ~  0,
  min_distance_scaled>0 ~ 1/min_distance_scaled
))

pu<-zero_one_scale(
  pu,
  cols = c("inv_dist", "inv_cost_pol")
)

# remove border patches:
pu$area<-pu$area/10^6
pu<-pu%>%filter(area >= 0.48713)


core_pa_for<-pa_optim(pu,target_lulc ="forest",lockin_col="core_prot_area_old",area_budget = gap[gap$lulc_name == "forest", ]$target_area/10^6)
core_pa_wat<-pa_optim(pu,target_lulc ="water",lockin_col="core_prot_area_old",area_budget = gap[gap$lulc_name == "water", ]$target_area/10^6)
core_pa_wetl<-pa_optim(pu,target_lulc ="wetland",lockin_col="core_prot_area_old",area_budget = gap[gap$lulc_name == "wetland", ]$target_area/10^6)
core_pa_agri<-pa_optim(pu,target_lulc ="agricultural",lockin_col="core_prot_area_old",area_budget = gap[gap$lulc_name == "agricultural", ]$target_area/10^6)
core_pa_global<-pa_optim(pu=pu,lockout = "built-up",lockin_col="core_prot_area_old",area_budget = gap[gap$lulc_name == "global_0.1_prot", ]$target_area/10^6)

PA_single_lulc<-rbind(core_pa_agri,core_pa_for,core_pa_wat,core_pa_wetl)

pu <- pu %>%
  left_join(PA_single_lulc %>%st_drop_geometry()%>% dplyr::select(ID, solution_1),
            by = "ID")
names(pu)[names(pu) == "solution_1"] <- "optim_PA_single_lulc"

pu <- pu %>%
  left_join(core_pa_global %>%st_drop_geometry()%>% dplyr::select(ID, solution_1),
            by = "ID")
names(pu)[names(pu) == "solution_1"] <- "optim_PA_global"


pu <- pu %>%
  mutate(core_pa_lulc = case_when(
    core_prot_area_old == TRUE            ~ "existing core pa",
    optim_PA_single_lulc == 1 & n_pa == 0 ~ "new core PA",
    optim_PA_single_lulc == 1 & n_pa >0 ~ "proposed upgrade existing PA",
    TRUE                       ~ "other"
  ))


pu <- pu %>%
  mutate(core_pa_global = case_when(
    core_prot_area_old == TRUE            ~ "existing core pa",
    optim_PA_global == 1 & n_pa == 0 ~ "new core PA",
    optim_PA_global == 1 & n_pa >0 ~ "proposed upgrade existing PA",
    TRUE                       ~ "other"
  ))

pu$core_pa_lulc <- factor(pu$core_pa_lulc,
                         levels = c(
                           "existing core pa",
                           "new core PA",
                           "proposed upgrade existing PA",
                           "other"
                         )
)

pu$core_pa_global <- factor(pu$core_pa_global,
                          levels = c(
                            "existing core pa",
                            "new core PA",
                            "proposed upgrade existing PA",
                            "other"
                          )
)

pu<-pu%>%mutate(is_core_prot_lulc = case_when(
  core_pa_lulc == "other" ~F,
  core_pa_lulc != "other" ~T
))

cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#00FF00",
  "proposed upgrade existing PA" = "#B8FFB8",
  "other" = NA
)

map_lulc <- ggplot(pu) +
  geom_sf(aes(fill = core_pa_lulc), color = NA) +
  scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(text = element_text(size = 20))

ggsave(paste0("WP4/2_output/02_optim/",siteID,"_pa_optim.png"), plot = map_lulc, width = 18, height = 10, dpi = 300)


stats_lulc_optim<-pu%>%mutate(core_pa = core_pa_lulc)%>%filter(core_pa!="other" )%>%group_by(lulc_name,core_pa)%>%
  summarise(ec_mean = mean(sampled_condition_scaled),
            ec_sd = sd(sampled_condition_scaled),
            km2 = sum(area))%>%st_drop_geometry() %>%mutate(scenario = "LULC")

stats_lulc_global<-pu%>%mutate(core_pa = core_pa_global)%>%filter(core_pa!="other" )%>%group_by(lulc_name,core_pa)%>%
  summarise(ec_mean = mean(sampled_condition_scaled),
            ec_sd = sd(sampled_condition_scaled),
            km2 = sum(area))%>%st_drop_geometry()%>%mutate(scenario = "global")


stats<-rbind(stats_lulc_optim,stats_lulc_global)

library(ggplot2)

ggplot(stats%>%filter(lulc_name != "built-up"),
       aes(x = scenario,
           y = km2,
           fill = core_pa)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap(~ lulc_name, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL,
    y = expression("Area (km"^2*")"),
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

pd <- position_dodge(width = 0.8)

ggplot(
  stats %>% filter(lulc_name != "built-up"),
  aes(x = scenario,
      y = ec_mean,
      fill = core_pa)
) +
  geom_col(
    position = pd,
    width = 0.7
  ) +
  geom_errorbar(
    aes(ymin = ec_mean - ec_sd,
        ymax = ec_mean + ec_sd),
    position = pd,
    width = 0.2,
    linewidth = 0.4
  ) +
  facet_wrap(~lulc_name, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL,
    y = "Mean ecosystem condition",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1))

# p_box <- ggplot(stats, aes(x = core_pa_lulc, y = sampled_habitat, fill = core_pa_lulc)) +
#   geom_boxplot() +
#   scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
#   theme_minimal() +
#   theme(legend.position = "none",text = element_text(size = 20))+
#   labs(
#     x = "",
#     y = "Ecosystem condition",
#     fill = ""
#   )

ggsave(paste0("WP4/2_output/02_optim/",siteID,"_ec_pa_optim.png"), plot = p_box, width = 18, height = 10, dpi = 300)

# kruskal.test(sampled_condition ~ core_prot, data = pu_stats)
# 
# pairwise.wilcox.test(pu_stats$sampled_condition, pu_stats$core_prot, p.adjust.method = "BH")

#save PU for OECM analysis
st_write(pu,paste0("WP4/2_output/02_optim/PA_optim_grid",siteID,".geojson"))
