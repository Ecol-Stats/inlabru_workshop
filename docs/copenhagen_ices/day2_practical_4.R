## ----child="practicals/spatial_data_types_areal.qmd"--------------------------

## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(CARBayesdata)

data(pollutionhealthdata)
data(GGHB.IZ)




## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(sf)
library(ggplot2)
library(scico)


## -----------------------------------------------------------------------------
resp_cases <- merge(GGHB.IZ, pollutionhealthdata, by = "IZ")


## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(dplyr)
resp_cases <- resp_cases %>% 
  mutate(SMR = observed/expected, .by = year )


## -----------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 6
#| fig-align: center
ggplot()+
  geom_sf(data=resp_cases,aes(fill=SMR))+
  facet_wrap(~year)+scale_fill_scico(palette = "roma")




## -----------------------------------------------------------------------------
#| message: false
#| warning: false

library(spdep)

W.nb <- poly2nb(GGHB.IZ,queen = TRUE)
W.nb




## -----------------------------------------------------------------------------
st_crs(GGHB.IZ)$units


## -----------------------------------------------------------------------------
W.nb250 <- poly2nb(GGHB.IZ,snap=250)
W.nb250




## -----------------------------------------------------------------------------
# subset the data
resp_cases_2011 <- resp_cases %>% filter(year ==2011)

# neighbors list 
nbw <- nb2listw(W.nb, style = "W")

# Global Moran's I
gmoran <- moran.test(resp_cases_2011$SMR, nbw,
                     alternative = "greater")
gmoran



## -----------------------------------------------------------------------------
lmoran <- localmoran(resp_cases_2011$SMR, nbw, alternative = "two.sided")


## -----------------------------------------------------------------------------

resp_cases_2011_m <- resp_cases_2011 %>% mutate(Zscores = lmoran[,"Z.Ii"])

resp_cases_2011_m <- resp_cases_2011_m %>% 
  mutate (SAC = case_when(
  Zscores > qnorm(0.975) ~ " M > 1" ,
  Zscores < -1*qnorm(0.975) ~ " M < 1",
  .default = "M = 0"),
  SAC=as.factor(SAC)
)



## -----------------------------------------------------------------------------
#| message: false
#| warning: false
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
library(mapview)
mapview(resp_cases_2011_m, 
        zcol = "SAC",
        layer.name = "SAC",
        col.regions=c("turquoise","orange","grey40"))



## ----child="practicals/spatial_data_types_geostats.qmd"-----------------------

## -----------------------------------------------------------------------------
#| message: false
#| warning: false

# For plotting
library(mapview)
library(ggplot2)
library(scico) # for colouring palettes

# Data manipulation
library(dplyr)




## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(sdmTMB)

pcod_df = sdmTMB::pcod 
qcs_grid = sdmTMB::qcs_grid



## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(sf)
pcod_sf =   st_as_sf(pcod_df, coords = c("lon","lat"), crs = 4326)


## -----------------------------------------------------------------------------
pcod_sf_proj <- st_transform(pcod_sf, crs = 32609)
st_crs(pcod_sf_proj)$units


## -----------------------------------------------------------------------------
pcod_sf_proj = st_transform(pcod_sf_proj,
                            gsub("units=m","units=km",
                                 st_crs(pcod_sf_proj)$proj4string)) 
st_crs(pcod_sf_proj)$units


## -----------------------------------------------------------------------------
pcod_sf = st_transform(pcod_sf,
                       crs = "+proj=utm +zone=9 +datum=WGS84 +no_defs +type=crs +units=km" )
st_crs(pcod_sf)$units



## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
#| message: false
#| warning: false
#| 
pcod_sf %>% 
  filter(year== 2017) %>%
  mutate(present = as.factor(present)) %>%
mapview(zcol = "present",
        layer.name = "Occupancy status of Pacific Cod in 2017")




## -----------------------------------------------------------------------------
#| message: false
#| warning: false

library(terra)
depth_r <- rast(qcs_grid, type = "xyz")
depth_r


## -----------------------------------------------------------------------------
crs(depth_r) <- crs(pcod_sf)


## -----------------------------------------------------------------------------
#| fig-width: 8 
#| fig-height: 8
#| fig-align: center  


library(tidyterra)

ggplot()+ 
  geom_spatraster(data=depth_r$depth)+
  geom_sf(data=pcod_sf,aes(color=factor(present))) +
  facet_wrap(~year)+
    scale_color_manual(name="Occupancy status for the Pacific Cod",
                     values = c("black","orange"),
                     labels= c("Absence","Presence"))+
  scale_fill_scico(name = "Depth",
                   palette = "nuuk",
                   na.value = "transparent" )





## -----------------------------------------------------------------------------
#| echo: false
#| fig-align: center
#| fig-width: 7
#| fig-height: 5

pcod_sf_subset <- pcod_sf %>% filter(year ==2017)

map <- ggplot(pcod_sf_subset) + geom_sf(aes(color = density), size = 2) +
scale_color_gradient(low = "blue", high = "orange") +
theme_bw()


point_1 <- st_buffer(pcod_sf_subset[1,],10)
point_2 <- st_buffer(pcod_sf_subset[12,],10)
# Extract coordinates and values for annotation
coords <- st_coordinates(pcod_sf_subset[c(1, 12), ])
values <- pcod_sf_subset$density[c(1, 12)]
# Offset the labels slightly to the right of the points
offset_x <- 10  # Adjust for horizontal shift
offset_y <- 10   # Adjust for vertical shift
# Create a DataFrame for text labels
labels_df <- data.frame(
  x = coords[,1]+ offset_x,
  y = coords[,2]+ offset_y,
  label = sprintf("Z(s) = %.2f\n(%.1f, %.1f)", values, coords[,1], coords[,2])
)

map+
  geom_sf(data=point_1,col=2,alpha=0.5)+
  geom_sf(data=point_2,col=4,alpha=0.5)+
  geom_label(data = labels_df, aes(x = x, y = y, label = label), 
             fill = "orange", # Soft background color
             color = "black", # Text color
             size = 3, 
             alpha=0.7,
             fontface = "bold",
             label.size = 0.3, # Mild border thickness
             label.r = unit(0.15, "lines"))




## -----------------------------------------------------------------------------
#| message: false
#| warning: false
#| fig-align: center
#| fig-width: 4
#| fig-height: 4


library(gstat)

pcod_sf_subset <- pcod_sf %>% filter(year ==2017)

vgm1 <- variogram(density~1, pcod_sf_subset)
plot(vgm1)



## -----------------------------------------------------------------------------
#| message: false
#| warning: false
#| fig-align: center
#| fig-width: 4
#| fig-height: 4

library(variosig)

varioEnv <- envelope(vgm1,
                     data = pcod_sf_subset,
                     locations = st_coordinates(pcod_sf_subset),
                     formula = density ~ 1,
                     nsim = 499)

envplot(varioEnv)




## ----child="practicals/spatial_data_types_points.qmd"-------------------------

## -----------------------------------------------------------------------------
#| message: false
#| warning: false

# For plotting
library(ggplot2)
library(scico) # for colouring palettes

# Data manipulation
library(dplyr)




## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(sf)
shp_SGC <-  st_read("datasets/SG_CairngormsNationalPark/SG_CairngormsNationalPark_2010.shp",quiet =T)



## -----------------------------------------------------------------------------
shp_SGC <- shp_SGC %>% st_transform(crs = 27700)
st_crs(shp_SGC)$units


## -----------------------------------------------------------------------------
shp_SGC <- st_transform(shp_SGC,gsub("units=m","units=km",st_crs(shp_SGC)$proj4string)) 
st_crs(shp_SGC)$units


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
ggplot()+
  geom_sf(data=shp_SGC)



## -----------------------------------------------------------------------------
ringlett <- read.csv("datasets/bnm_ringlett.csv")
head(ringlett)


## -----------------------------------------------------------------------------
ringlett_sf <- ringlett %>% st_as_sf(coords = c("x","y"),crs = "+proj=longlat +datum=WGS84") 




## -----------------------------------------------------------------------------
ringlett_CNP <- ringlett_sf[shp_SGC,] # crop to mainland


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
ggplot()+
  geom_sf(data=shp_SGC)+
  geom_sf(data=ringlett_CNP)


## -----------------------------------------------------------------------------
#| message: false
#| warning: false
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
library(terra)
elevation_r <- rast("datasets/Scotland_elev.tiff")
crs(elevation_r) = crs(shp_SGC)
plot(elevation_r)


## -----------------------------------------------------------------------------
elevation_r <- elevation_r %>% scale()


## -----------------------------------------------------------------------------
#| fig-width: 6
#| fig-height: 4
#| fig-align: center

elev_CNP <- terra::crop(elevation_r,shp_SGC,mask=T)
plot(elev_CNP)



