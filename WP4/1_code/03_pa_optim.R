# test prioritizr
library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(gdistance)
source("WP4/1_code/wp4_functions_utils.R")

main_dir<-"P:/312204_pareus/"
siteID<-"FRL04"
gap<-read.csv(paste0("WP4/2_output/01_PA_analysis/",siteID,"_gap_analysis.csv"))

stud_area<-read_sf(paste0(main_dir,"WP2/T2.2/PGIS_ES_mapping/",siteID,"/raw_data_backup/study_site.gpkg"))
stud_area<-stud_area%>%filter(siteID=="FRL04")

pu<-st_read(paste0("WP4/2_output/02_optim/",siteID,"_input_final_grid.json"))
pu$ID<-c(1:nrow(pu))
pu$area<-as.numeric(st_area(pu))
pu$inv_cost_pol<-1/pu$sampled_cost_pol

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
    core_prot_area_old == TRUE   ~ "existing core pa",
    optim_PA_single_lulc == 1 & n_pa == 0 ~ "new core PA",
    optim_PA_single_lulc == 1 & n_pa >0 ~ "proposed upgrade existing PA",
    core_prot_area_old == FALSE & is.na(optim_PA_single_lulc) & n_pa == 0 ~ "not protected",
    core_prot_area_old == FALSE & is.na(optim_PA_single_lulc) & n_pa > 0 ~ "other protected areas"
  ))


pu <- pu %>%
  mutate(core_pa_global = case_when(
    core_prot_area_old == TRUE  ~ "existing core pa",
    optim_PA_global == 1 & n_pa == 0 ~ "new core PA",
    optim_PA_global == 1 & n_pa >0 ~ "proposed upgrade existing PA",
    core_prot_area_old == FALSE & is.na(optim_PA_global)  & n_pa == 0 ~ "not protected",
    core_prot_area_old == FALSE & is.na(optim_PA_global)  & n_pa > 0 ~ "other protected areas"
  ))

pu$core_pa_lulc <- factor(pu$core_pa_lulc,
                         levels = c(
                           "existing core pa",
                           "new core PA",
                           "proposed upgrade existing PA",
                           "not protected",
                           "other protected areas"
                         )
)

pu$core_pa_global <- factor(pu$core_pa_global,
                            levels = c(
                              "existing core pa",
                              "new core PA",
                              "proposed upgrade existing PA",
                              "not protected",
                              "other protected areas"
                            )
)

pu<-pu%>%mutate(is_core_prot_lulc = case_when(
  core_pa_lulc == "other" ~F,
  core_pa_lulc != "other" ~T
))

pu<-pu%>%mutate(is_core_prot_global = case_when(
  core_pa_global == "other" ~F,
  core_pa_global != "other" ~T
))

cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#bf8bff",
  "proposed upgrade existing PA" = "#e5d0ff",
  "other" = NA
)

map_lulc <- ggplot(pu) +
  geom_sf(aes(fill = core_pa_lulc), color = NA) +
  scale_fill_manual(values = cols, name = NULL, na.translate = FALSE) +
  geom_sf(data = stud_area, fill = NA, color = "black") +
  theme_minimal()+
  theme(
    legend.position = "top",
    text = element_text(size = 15)
  )+ggtitle("PA single LULC")

ggsave(paste0("WP4/2_output/02_optim/",siteID,"_pa_global_optim.png"), plot = map_lulc, width = 18, height = 10, dpi = 300)


stats_lulc_optim<-pu%>%mutate(core_pa = core_pa_lulc)%>%group_by(lulc_name,core_pa)%>%
  summarise(ec_mean = mean(sampled_condition_scaled),
            ec_sd = sd(sampled_condition_scaled),
            km2 = sum(area),
            reg_mean = mean(sampled_reg_scaled),
            reg_sd = sd(sampled_reg_scaled),
            pol_cost_mean = mean(sampled_cost_pol),
            cost_sd = sd(sampled_cost_pol),
            dist_mean = mean(min_distance),
            dist_sd = sd(min_distance))%>%st_drop_geometry() %>%mutate(scenario = "LULC")

stats_global_optim<-pu%>%mutate(core_pa = core_pa_global)%>%group_by(lulc_name,core_pa)%>%
  summarise(ec_mean = mean(sampled_condition_scaled),
            ec_sd = sd(sampled_condition_scaled),
            km2 = sum(area),
            reg_mean = mean(sampled_reg_scaled),
            reg_sd = sd(sampled_reg_scaled),
            pol_cost_mean = mean(sampled_cost_pol),
            cost_sd = sd(sampled_cost_pol),
            dist_mean = mean(min_distance),
            dist_sd = sd(min_distance))%>%st_drop_geometry()%>%mutate(scenario = "GLOBAL")


stats<-rbind(stats_lulc_optim,stats_global_optim)

#area statistics
cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#bf8bff",
  "proposed upgrade existing PA" = "#e5d0ff",
  "other protected areas"        = "#ADD8E6",
  "not protected"                = "#1A1A1A"
)
ggplot(stats%>%filter(core_pa != "other protected areas" & !is.na(lulc_name) & lulc_name != "built-up" & core_pa != "not protected" ),
       aes(x = scenario,
           y = km2,
           fill = core_pa)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap(~ lulc_name, scales = "free_y") +
  scale_fill_manual(values = cols)+
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


st_write(pu,paste0("WP4/2_output/02_optim/PA_optim_grid",siteID,".geojson"))

cols <- c(
  "existing core pa" = "#00A300",
  "new core PA" = "#bf8bff",
  "proposed upgrade existing PA" = "#e5d0ff",
  "other protected areas"        = "#ADD8E6",
  "not protected"                = "#1A1A1A"
)

ltypes <- c(
  "existing core pa"             = "dashed",
  "proposed upgrade existing PA" = "solid",
  "new core PA"                  = "solid",
  "other protected areas"        = "solid",
  "not protected"                = "solid"
)

stats %>%
  filter(!is.na(lulc_name)) %>%
  mutate(
    lulc_name = factor(lulc_name,
                       levels = unique(lulc_name))
  ) %>%
  ggplot(
    aes(
      x = lulc_name,
      y = dist_mean/1000,
      group = core_pa,
      colour = core_pa,
      linetype = core_pa
    )
  ) +
  geom_path(linewidth = 1.3) +
  geom_point(size = 4) +
  scale_colour_manual(values = cols)+
  scale_linetype_manual(values = ltypes) +
  coord_polar() +
  facet_wrap(~scenario) +
  # scale_y_continuous(
  #   limits = c(0, 1),
  #   breaks = seq(0, 1, 0.2)
  # ) +
  labs(
    x = NULL,
    y = "Distance to existing core PA [km] <-- min",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  )+
  guides(
    colour = guide_legend(override.aes = list(
      linetype = unname(ltypes),
      linewidth = 1.2
    )),
    linetype = "none"
  )





