library(giscoR)
library(sf)
library(dplyr)

# Download French NUTS-3 polygons
fr_dep <- gisco_get_nuts(
  year = "2021",
  country = "FR",
  nuts_level = "3",
  epsg = "4326",
  resolution = "01"
)

# Select departments
dep <- fr_dep %>%
  filter(NUTS_ID %in% c(
    "FRL01",  # Alpes-de-Haute-Provence
    "FRL02",  # Hautes-Alpes
    #"FRL03",  # Alpes-Maritimes
    "FRL04",  # Bouches-du-Rhône
    #"FRL05", #Var
    "FRK23",  #Drome
    "FRL06" # Vaucluse
  ))


# communes <- gisco_get_communes(
#   year = 2016,
#   country = "FR",
#   epsg = 3035,
#   spatialtype = "RG"
# )

# Make sure your geometry uses the same CRS
# geom <- st_transform(geom, 3035)



cases<-st_read("data/shared/pareus_sites.gpkg")
cases_sel<-cases%>%filter(siteID %in% c("FRA_BAR2","FRL04"))
# FRA_BAR2<-cases%>%filter(siteID %in% c("FRA_BAR2"))
# Select communes that intersect the geometry
# communes <- st_transform(communes, st_crs(FRA_BAR2))
# communes_sel <- communes[
#   lengths(st_intersects(communes, FRA_BAR2)) > 0,
# ]
# plot(st_geometry(communes_sel))
plot(st_geometry(dep), border = "blue")

plot(st_geometry(cases),border = "red",add=T, lwd = 2)
text(
  st_coordinates(st_centroid(dep)),
  labels = dep$NAME_LATN,
  cex = 0.8,
  col = "black"
)


dep_dissolved <- st_union(dep)
st_write(dep_dissolved,"data/shared/FRA.gpkg")
st_write(dep_dissolved,"data/shared/FRA.shp")



### NOR
NOR <- gisco_get_communes(
  year = 2016,
  country = "NOR",
  epsg = 3035,
  spatialtype = "RG"
)
cases_sel<-cases%>%filter(siteID %in% c("TRD"))
NOR_sel<-NOR%>%filter(COMM_ID %in% c("NO1663","NO1657","NO1662","NO1653","NO1601"))
selb<-NOR%>%filter(COMM_ID %in% c("NO1664"))

NOR_sel<-st_transform(NOR_sel,st_crs(cases))
selb<-st_transform(selb,st_crs(cases))

selb_valid <- st_make_valid(selb)
cases_sel <- st_make_valid(cases_sel)

# Run the intersection on the cleaned data
sf_use_s2(FALSE)
a <- st_intersection(selb_valid, cases_sel)
a<-a[,c(1:10)]
plot(st_geometry(NOR_sel))
text(
  st_coordinates(st_centroid(NOR_sel)),
  labels = NOR_sel$NAME_LATN,
  cex = 0.8,
  col = "black"
)

plot(st_geometry(cases_sel),border = "red",add=T, lwd = 2)
plot(st_geometry(a),border = "blue",add=T, lwd = 2)

NOR_sel<-rbind(a,NOR_sel)

NOR_sel_diss <- st_union(NOR_sel)
st_write(NOR_sel_diss,"data/shared/NOR.gpkg")
st_write(NOR_sel_diss,"data/shared/NOR.shp")
