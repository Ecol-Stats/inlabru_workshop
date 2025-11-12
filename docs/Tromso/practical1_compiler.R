## ----child="practicals/distance_sampling.qmd"---------------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     
library(sf)
# load some libraries to generate nice map plots
library(scico)
library(mapview)


## -----------------------------------------------------------------------------
mexdolphin <- mexdolphin_sf
mexdolphin$depth <- mexdolphin$depth %>% mutate(depth=scale(depth)%>%c())
mapviewOptions(basemaps = c( "OpenStreetMap.DE"))

mapview(mexdolphin$points,zcol="size")+
  mapview(mexdolphin$samplers)+
 mapview(mexdolphin$ppoly )





## -----------------------------------------------------------------------------
boundary0 = fm_nonconvex_hull(mexdolphin$points,convex = -0.1)

mesh_0 = fm_mesh_2d(boundary = boundary0,
                          max.edge = c(30, 150), # The largest allowed triangle edge length.
                          cutoff = 15,
                          crs = fm_crs(mexdolphin$points))
ggplot() + gg(mesh_0)


## -----------------------------------------------------------------------------
mesh_1 = fm_mesh_2d(boundary = mexdolphin$ppoly,
                    max.edge = c(30, 150),
                    cutoff = 15,
                    crs = fm_crs(mexdolphin$points))
ggplot() + gg(mesh_1)




## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(terra)
library(stars)
# Convert sf to stars raster
stars_raster <- st_rasterize(mexdolphin$depth[, "depth"])
# Convert stars to terra raster if needed
terra_raster <- rast(stars_raster)


## -----------------------------------------------------------------------------
# Extend raster ext by 5 % of the original raster
re <- extend(terra_raster, ext(terra_raster)*2.1)
# Convert to an sf spatial object
re_df <- re %>% stars::st_as_stars() %>%  st_as_sf(na.rm=F)
# fill in missing values using the original raster 
re_df$depth <- bru_fill_missing(terra_raster,re_df,re_df$depth)
# store the projectes values as a raster
depth_rast_p <- stars::st_rasterize(re_df) %>% rast()


## -----------------------------------------------------------------------------
#| echo: false
#| fig-align: center
#| message: false
#| warning: false
library(tidyterra)

ggplot()+geom_sf(data=mexdolphin$depth,aes(color=depth))+scale_color_scico(palette = "roma")+ggtitle("Depth sf object") +gg(mexdolphin$mesh) + 
ggplot()+tidyterra::geom_spatraster(data=depth_rast_p)+scale_fill_scico(palette = "roma")+ggtitle("Depth raster object")+gg(mexdolphin$mesh)



## -----------------------------------------------------------------------------
spde_model <- inla.spde2.pcmatern(mexdolphin$mesh,
  prior.sigma = c(2, 0.01),
  prior.range = c(50, 0.01)
)




## -----------------------------------------------------------------------------
hn <- function(distance, sigma) {
  exp(-0.5 * (distance / sigma)^2)
}


## -----------------------------------------------------------------------------
cmp <- ~ space(main = geometry, model = spde_model) +
  sigma(1,
    prec.linear = 1,
    marginal = bm_marginal(qexp, pexp, dexp, rate = 1 / 8)
  ) +
  depth(depth_rast_p$depth,model="linear")+
  Intercept(1)


## -----------------------------------------------------------------------------
#| echo: true
#| eval: false
#| 
# eta <- ... + log(2)




## -----------------------------------------------------------------------------
# build integration scheme
distance_domain <-  fm_mesh_1d(seq(0, 8,
                              length.out = 30))
ips = fm_int(list(geometry = mexdolphin$mesh,
                  distance = distance_domain),
             samplers = mexdolphin$samplers)


## -----------------------------------------------------------------------------
lik = bru_obs("cp",
              formula = eta,
              data = mexdolphin$points,
              ips = ips)


## -----------------------------------------------------------------------------
fit = bru(cmp, lik)






## -----------------------------------------------------------------------------
plot( spde.posterior(fit, "space", what = "range")) +
plot( spde.posterior(fit, "space", what = "log.variance"))  


## -----------------------------------------------------------------------------
pxl <- fm_pixels(mexdolphin$mesh, dims = c(200, 100), mask = mexdolphin$ppoly)


## -----------------------------------------------------------------------------
pr.int = predict(fit, pxl, ~data.frame(spatial = space,
                                      lambda = exp(Intercept + depth + space)))


## -----------------------------------------------------------------------------
ggplot() + geom_sf(data = pr.int$spatial,aes(color = mean)) + scale_color_scico() + ggtitle("Posterior mean")

ggplot() + geom_sf(data = pr.int$spatial,aes(color = sd)) + scale_color_scico() + ggtitle("Posterior sd")




## -----------------------------------------------------------------------------
distdf <- data.frame(distance = seq(0, 8, length.out = 100))
dfun <- predict(fit, distdf, ~ hn(distance, sigma))
plot(dfun)


## -----------------------------------------------------------------------------
predpts <- fm_int(mexdolphin$mesh, mexdolphin$ppoly)
Lambda <- predict(fit, predpts, ~ sum(weight * exp(space + Intercept)))
Lambda


## -----------------------------------------------------------------------------
Ns <- seq(50, 450, by = 1)
Nest <- predict(fit, predpts,
  ~ data.frame(
    N = Ns,
    density = dpois(
      Ns,
      lambda = sum(weight * exp(space + Intercept))
    )
  ),
  n.samples = 2000
)


## -----------------------------------------------------------------------------
Nest <- dplyr::bind_rows(
  cbind(Nest, Method = "Posterior"),
  data.frame(
    N = Nest$N,
    mean = dpois(Nest$N, lambda = Lambda$mean),
    mean.mc_std_err = 0,
    Method = "Plugin"
  )
)


## -----------------------------------------------------------------------------
ggplot(data = Nest) +
  geom_line(aes(x = N, y = mean, colour = Method)) +
  geom_ribbon(
    aes(
      x = N,
      ymin = mean - 2 * mean.mc_std_err,
      ymax = mean + 2 * mean.mc_std_err,
      fill = Method,
    ),
    alpha = 0.2
  ) +
  geom_line(aes(x = N, y = mean, colour = Method)) +
  ylab("Probability mass function")


## -----------------------------------------------------------------------------
bc <- bincount(
  result = fit,
  observations = mexdolphin$points$distance,
  breaks = seq(0, max(mexdolphin$points$distance), length.out = 9),
  predictor = distance ~ hn(distance, sigma)
)
attributes(bc)$ggp


## -----------------------------------------------------------------------------
hr <- function(distance, sigma) {
  1 - exp(-(distance / sigma)^-1)
}



## ----child="practicals/distance_sampling_spat.qmd"----------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     
library(sf)
# load some libraries to generate nice map plots
library(scico)
library(mapview)


## -----------------------------------------------------------------------------

mrsea <- inlabru::mrsea

ggplot() +
  geom_fm(data = mrsea$mesh) +
  gg(mrsea$boundary) +
  gg(mrsea$samplers) +
  gg(mrsea$points, size = 0.5) +
  facet_wrap(~season) +
  ggtitle("MRSea observation seasons")


## -----------------------------------------------------------------------------
matern <- inla.spde2.pcmatern(mrsea$mesh,
  prior.sigma = c(0.1, 0.01),
  prior.range = c(10, 0.01)
)


## -----------------------------------------------------------------------------
cmp <- ~ Intercept(1) + 
  space_time(
    geometry,
    model = matern,
    group = season,
    ngroup = 4,
    control.group = list(model="iid")
  )


## -----------------------------------------------------------------------------
#| eval: false
# eta <- ... + ... ~ ... + ...




## -----------------------------------------------------------------------------
ips <- fm_int(
  domain = list(geometry = mrsea$mesh, season = 1:4),
  samplers = mrsea$samplers
)


## -----------------------------------------------------------------------------
#| eval: false
# lik = bru_obs(formula = ...,
#     family = ...,
#     data = ...,
#     ips  = ...)
# fit = bru(..., ...)
# 




## -----------------------------------------------------------------------------
ppxl <- fm_pixels(mrsea$mesh, mask = mrsea$boundary, format = "sf")
ppxl_all <- fm_cprod(ppxl, data.frame(season = seq_len(4)))

lambda1 <- predict(
  fit,
  ppxl_all,
  ~ data.frame(season = season, lambda = exp(space_time + Intercept))
)

pl1 <- ggplot() +
  gg(lambda1, geom = "tile", aes(fill = q0.5)) +
  gg(mrsea$points, size = 0.3) +
  facet_wrap(~season) +
  coord_sf()
pl1


