#### case study plot
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# 1. Install and load the modern package
#install.packages("giscoR")
library(giscoR)

# Download national boundary of France
sites <- gisco_get_nuts(year = "2024", nuts_id = c("FRL0","SK021","NO06"))
# 
# main_dir<-"P:/312204_pareus/"
# 
# sites<-sf::st_read(paste0(main_dir,"pareus_repository/shared/pareus_sites.gpkg"))

sites<-st_make_valid(sites)
# Label positions (centroids)
sf::sf_use_s2(FALSE)

labels <- sites %>%
  st_centroid() %>%
  mutate(
    x = c(11.2,  6, 18),
    y = c(65.5, 45.7, 49.5)
  )

sf::sf_use_s2(TRUE)

# Europe
europe <- ne_countries(
  continent = "Europe",
  scale = "medium",
  returnclass = "sf"
)

# sites<-sites%>%filter(siteID %in% c("SK021","TRD","FRA"))
# sites<-st_collection_extract(sites, "POLYGON")

ggplot() +
  # Sea
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "#dceef8"),
    panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_blank()
  ) +
  
  # Countries
  geom_sf(
    data = europe,
    fill = "grey96",
    colour = "grey60",
    linewidth = 0.25
  ) +
  
  # Study areas
  geom_sf(
    data = sites,
    fill = "#D55E00",
    colour = "black",
    linewidth = 0.5,
    alpha = 0.8
  ) +
  
  # Labels
  geom_text(
    data = labels,
    aes(x, y, label = CNTR_CODE),
    fontface = "bold",
    size = 3.5
  ) +
  
  coord_sf(
    xlim = c(-12, 28),
    ylim = c(42, 71),
    expand = FALSE,
    crs = 4326,
    label_graticule = "SW"
  ) +
  
  scale_x_continuous(breaks = seq(-10, 30, 5)) +
  scale_y_continuous(breaks = seq(40, 70, 5))

