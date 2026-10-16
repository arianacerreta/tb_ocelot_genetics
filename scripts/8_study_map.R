####study map

##libraries
library("ggplot2")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("terra")
library("tidyterra")
library("raster")
library("tidyverse")
library("prismatic")
library("ggmap")
library("patchwork")
library("cowplot")
library("ggspatial")

#layers
#countries
world<- ne_countries(scale = "medium", returnclass = "sf")

background<- ne_download(scale = 50, type ="HYP_50M_SR_W", category = "raster")
#states
states<-ne_download(scale = 50L, type = "states", category = "cultural")

usstate<-st_read("./data/spatial_inputs/tl_2026_us_state/tl_2026_us_state.shp")
#refuge properties
NWR<-st_read("./data/spatial_inputs/National_Wildlife_Refuge_System_Boundaries/FWSBoundaries.shp",
             type = 1)
#ranch
elsauz<-st_read("./data/spatial_inputs/ElSauz_boundary/ElSauz_boundary/ElSauz_boundary.shp")

#IUCN ocelot distribution
ocelot_dist<-st_read("./data/spatial_inputs/IUCN_ocelot_dist/data_0.shp")

#historical
historical<-st_read("./data/spatial_inputs/HistoricRange/HistoricRange/Haines_HistoricRange.shp")
#subset layers
americas<-subset(world, continent %in% c("North America", "South America") |
                   geounit == "French Guiana")

france<-subset(world, subunit %in% c("France"))
fr_tracts<-st_cast(france, "POLYGON")
frenchguiana<-st_geometry(fr_tracts[c(10),])

TX<-subset(states, iso_3166_2 %in% c("US-TX"))

LANWR<-subset(NWR, ORGCODE %in% c("21553")) #USFWS requested using entire LANWR boundary
LRGVNWR<-subset(NWR, ORGCODE %in% c("21552"))

states_tract<-st_cast(usstate, "POLYGON")
TX_tl<-st_geometry(states_tract[c(79),])

#pull just Yturria
tracts<-st_cast(LRGVNWR, "POLYGON")
yturria<-st_geometry(tracts[c(170:171),])

#crop background
extent<-ext(-120, -34,-31,39)
background_cropped<-crop(background, extent, ext = FALSE)

#check CRS
st_crs(americas) #epsg 4326; WGS84
st_crs(ocelot_dist) #espg 4326; WGS84
st_crs(TX) #espg 4326; WGS84
st_crs(LANWR) #espg 4269; NAD83 
st_crs(yturria) #espg 4269; NAD83 
st_crs(elsauz)# EPSG 26914
st_crs(background_cropped) #espg 4326
st_crs(TX_tl)#espg 4269; NAD83 
#distribution map

#project texas in NAD33 for plotting
TX_NAD83<-st_transform(TX, crs = st_crs(4269))
americas_NAD83<-st_transform(americas, crs = st_crs(4269))
elsauz_NAD83<-st_transform(elsauz, crs = st_crs(4269))
yturria_NAD83<-st_transform(yturria, crs = st_crs(4269))
LANWR_NAD83<-st_transform(LANWR,crs = st_crs(4269))

crop_extent<-extent(-97.8, -97.1, 25.9,26.7)

LANWR_WGS84<-st_transform(LANWR,crs = st_crs(4326))



#landsat?

filenames<-list.files("./data/spatial_inputs/landsat/", pattern = "*.TIF", full.names = T)
combined_landsat<-rast(filenames)
crs(combined_landsat) #PROJCRS[\"WGS 84 / UTM zone 14N\ (epsg:32614)

ggplot()+
  geom_spatraster(data=combined_landsat)+
  labs(fill="")+
  facet_wrap(~lyr)+
  scale_fill_hypso_c("colombia",na.value = NA)+
 # geom_sf(data = TX_tl, fill = NA, color = "#1d2e28")+
  theme_classic()

crop_extent<-extent(600000, 685000, 2875000,2950000)
ext(combined_landsat)
cropped_landsat<-crop(combined_landsat, crop_extent)

#SpatExtent : 543585, 770715, 2759385, 2991315 (xmin, xmax, ymin, ymax)
#projected
TX_32614<-st_transform(TX_tl, crs = st_crs(32614))
americas_32614<-st_transform(americas, crs = st_crs(32614))
elsauz_32614<-st_transform(elsauz, crs = st_crs(32614))
yturria_32614<-st_transform(yturria, crs = st_crs(32614)) #xmin: 634511.9 ymin: 2940689 xmax: 637805.7 ymax: 2942388
LANWR_32614<-st_transform(LANWR,crs = st_crs(32614)) # xmin: 657766.3 ymin: 2878694 xmax: 680764.3 ymax: 2939301


ranch<-st_union(elsauz_32614, yturria_32614)
#"ranch" = "#D55E00", "refuge" = "#0072B2"

main_map<-ggplot()+
  geom_spatraster_rgb(data=cropped_landsat,
                      r=4,g=3,b=2,
                      max_col_value = 30000,
                      interpolate=T)+
 # geom_sf(data = TX_32614, fill = NA, color = "#1d2e28")+
  geom_sf(data = ranch, aes(fill = "Ranch"), color = "#bc5300", alpha = 0.4, size = 0.7)+ #ranch
  geom_sf(data = LANWR_32614, aes(fill = "Refuge"), color = "#00517f", alpha = 0.6,  size = 0.7)+
  coord_sf(xlim = c(600000, 685000), ylim = c(2875000,2950000), expand = FALSE)+
  scale_fill_manual(values = c("Ranch" = "#D55E00", "Refuge" = "#0072B2"))+
  theme_void()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.box.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        legend.position = c(0,0.81),
        legend.margin = margin(t=0.25, b= 0.25, l=0.25, r =0.25, unit = "cm"),
        legend.justification = c("left", "top"))+
  annotation_scale(data= cropped_landsat, location = "br")
  

#make extent of main map into polygon; project to #epsg 4326; WGS84
inset_extent<-project(cropped_landsat, y= "epsg:4326")

inset<-ggplot()+
  geom_spatraster_rgb(data = background_cropped, maxcell = 5e6)+
  geom_sf(data = americas, fill = "gray", color = "gray50")+
  geom_sf(data = frenchguiana, fill = "gray", color = "gray50")+
  geom_sf(data = TX, fill = NA, color = "gray50")+
  geom_sf(data = ocelot_dist, fill = "#742e8f", color = "#742e8f", alpha = 0.5)+
  coord_sf(xlim = c(-111.2, -34.5), ylim = c(-30.6,36.9), expand = FALSE)+
  theme_void()+
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1, "cm"),
                           width = unit(1, "cm"))+
  geom_rect(data = data.frame(xmin = st_bbox(inset_extent)[[1]],
                              ymin = st_bbox(inset_extent)[[2]],
                              xmax = st_bbox(inset_extent)[[3]],
                              ymax = st_bbox(inset_extent)[[4]]),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA, 
    colour = "black",
    linewidth = 0.5, inherit.aes = FALSE
  )

pretty_inset<-ggdraw(inset)+
  draw_line(x=c(0,0,1,1,0), y=c(0.005,0.995,0.995,0.005,0.005), color = "black", size = 0.7)+
  draw_label("Current Ocelot Distribution", x=0.7, y=0.35, color = "black", size = 11)

ggdraw(main_map)+
  draw_plot(pretty_inset, x = 0, y = 0.004, width = 0.67, height = 0.67)
  
