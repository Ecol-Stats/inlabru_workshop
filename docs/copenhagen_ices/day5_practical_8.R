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
  Intercept(1)


## -----------------------------------------------------------------------------
eta <- geometry + distance ~ space +
  log(hn(distance, sigma)) +
  Intercept + log(2) 


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
                                      lambda = exp(Intercept + space)))


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



## ----child="practicals/zero_inflated_ex.qmd"----------------------------------

## -----------------------------------------------------------------------------
#| echo: false
#| include: false
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE,
  dev = "png",
  dev.args = list(type = "cairo-png"),
  fig.width = 7,
  fig.height = 5
)

# load webexercises library for tasks and questions (just for a preview - the practical compiler should take care of this when compiling multiple excercises)
library(webexercises)


## ----setup--------------------------------------------------------------------
#| eval: true
#| echo: true
#| include: true
#| message: false
#| warning: false
#| label: setup

library(dplyr)
library(ggplot2)
library(inlabru)
library(INLA)
library(terra)
library(sf)
library(scico)
library(magrittr)
library(patchwork)
library(tidyterra)


# We want to obtain CPO data from the estimations
bru_options_set(control.compute = list(dic = TRUE,
                                       waic = TRUE,
                                       mlik = TRUE,
                                       cpo = TRUE))


## -----------------------------------------------------------------------------
#| fig-cap: "Location of gorilla nests"
#| out-width: "80%"
#| fig-align: 'center'
#| label: nests_loc
gorillas_sf <- inlabru::gorillas_sf
nests <- gorillas_sf$nests
boundary <- gorillas_sf$boundary

ggplot() + geom_sf(data = nests) +
  geom_sf(data = boundary, alpha = 0)



## -----------------------------------------------------------------------------
gcov = gorillas_sf_gcov()
elev_cov <- gcov$elevation
dist_cov <-  gcov$waterdist


## -----------------------------------------------------------------------------
#| fig-cap: "Covariates"
#| out-width: "80%"
#| fig-align: 'center'
#| echo: false
#| 
theme_map = theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

p1 = ggplot() + geom_spatraster(data = elev_cov) + ggtitle("Elevation") +
  scale_fill_scico(direction = -1) + geom_sf(data = boundary, alpha = 0) + theme_map
p2 = ggplot() + geom_spatraster(data = dist_cov) + ggtitle("Distance to water") +
  scale_fill_scico(direction = -1) + geom_sf(data = boundary, alpha = 0) +
  theme_map
  
p1 + p2


## -----------------------------------------------------------------------------
#| fig-cap: "Counts of gorilla nests"
#| out-width: "80%"
#| fig-align: 'center'
#| 
# Rasterize data
counts_rstr <-
  terra::rasterize(vect(nests), gcov, fun = sum, background = 0) %>%
  terra::aggregate(fact = 5, fun = sum) %>%
  mask(vect(sf::st_geometry(boundary)))
plot(counts_rstr)
# compute cell area
counts_rstr <- counts_rstr %>%
  cellSize(unit = "km") %>%
  c(counts_rstr)


## -----------------------------------------------------------------------------
#| echo: true
counts_df <- crds(counts_rstr, df = TRUE, na.rm = TRUE) %>%
  bind_cols(values(counts_rstr, mat = TRUE, na.rm = TRUE)) %>%
  rename(count = sum) %>%
  mutate(present = (count > 0) * 1L) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(nests))


## -----------------------------------------------------------------------------
elev_cov1 <- elev_cov %>% 
  terra::aggregate(fact = 5, fun = mean) %>% scale()
dist_cov1 <- dist_cov %>% 
  terra::aggregate(fact = 5, fun = mean) %>% scale()


## -----------------------------------------------------------------------------
#| label: fig-covariate-raster
#| fig-cap: "Covariates "
#| out-width: "80%"
#| fig-align: 'center'
#| echo: false

p1 = ggplot() + geom_spatraster(data = elev_cov1) + ggtitle("Elevation") +
  scale_fill_scico(direction = -1) + geom_sf(data = boundary, alpha = 0) + theme_map
p2 = ggplot() + geom_spatraster(data = dist_cov1) + ggtitle("Distance to water") +
  scale_fill_scico(direction = -1) + geom_sf(data = boundary, alpha = 0) +
  theme_map
  
p1 + p2



## -----------------------------------------------------------------------------

mesh <- fm_mesh_2d(
  loc = st_as_sfc(counts_df),
  max.edge = c(0.5, 1),
  crs = st_crs(counts_df)
)

matern <- inla.spde2.pcmatern(mesh,
  prior.sigma = c(1, 0.01),
  prior.range = c(5, 0.01)
)


## -----------------------------------------------------------------------------
#| fig-cap: "Mesh over the count locations"
#| fig-align: 'center'
#| out-width: "80%"
#| echo: false
ggplot() +
  geom_fm(data = mesh) +
  geom_sf(
    data = counts_df[counts_df$count > 0, ],
    aes(color = count),
    size = 1,
    pch = 4
  ) +
  theme_minimal() + theme_map


## -----------------------------------------------------------------------------
#| echo: true
#| eval: false

# cmp = ~ Intercept(1) + elevation(...) + distance(...) + space(...)
# 
# lik = bru_obs(...,
#     E = area)
# 
# fit_zip <- bru(cmp, lik)






## -----------------------------------------------------------------------------
#| fig-cap: "Estimated $\\lambda$ (left) and expected counts (right) with zero inflated model"
#| fig-align: 'center'
#| out-width: "80%"
#| echo: false
p1 = ggplot() + gg(pred_zip$lambda, aes(fill = mean), geom="tile") + theme_map + scale_fill_scico(direction = -1) + 
  ggtitle(expression("Posterior mean of " ~lambda)) + theme(legend.position = "bottom")
p2 = ggplot() + gg(pred_zip$expect, aes(fill = mean), geom="tile") + theme_map + scale_fill_scico(direction = -1)+ 
  ggtitle(expression("Posterior mean of Expected counts"))+ theme(legend.position = "bottom")
p1 + p2







## -----------------------------------------------------------------------------
#| fig-cap: "Estimated $\\lambda$ (left) and expected counts (right) with hurdle model"
#| fig-align: 'center'
#| out-width: "80%"
#| echo: false
p1 = ggplot() + gg(pred_zap$lambda, aes(fill = mean), geom="tile") + theme_map + scale_fill_scico(direction = -1) + 
  ggtitle(expression("Posterior mean of " ~lambda)) + theme(legend.position = "bottom")
p2 = ggplot() + gg(pred_zap$expect, aes(fill = mean), geom="tile") + theme_map + scale_fill_scico(direction = -1)+ 
  ggtitle(expression("Posterior mean of Expected counts"))+ theme(legend.position = "bottom")
p1 + p2



## -----------------------------------------------------------------------------
#| echo: true
#| eval: false

# # define components
# cmp <- ~
#   Intercept_count(1) +
#     elev_count(elev_cov1, model = "linear") +
#     dist_count(dist_cov1, model = "linear") +
#     space_count(geometry, model = matern) +
#     Intercept_presence(1) +
#     elev_presence(elev_cov1, model = "linear") +
#     dist_presence(dist_cov1, model = "linear") +
#     space_presence(geometry, model = matern)
# 
# # positive count model
# pos_count_obs <- bru_obs(formula = ...,
#       family = ...,
#       data = counts_df[counts_df$present > 0, ],
#       E = area)
# 
# # presence model
# presence_obs <- bru_obs(formula ...,
#   family = ...,
#   data = counts_df,
# )
# 
# # fit the model
# fit_zap2 <- bru(...)




## ----fig-fit-zap--------------------------------------------------------------
#| eval: true
#| echo: true

cmp <- ~
  Intercept_count(1) +
    elev_count(elev_cov1, model = "linear") +
    dist_count(dist_cov1, model = "linear") +
    Intercept_presence(1) +
    elev_presence(elev_cov1, model = "linear") +
    dist_presence(dist_cov1, model = "linear") +
    space(geometry, model = matern) +
  space_copy(geometry, copy = "space", fixed = FALSE)


pos_count_obs <- bru_obs(formula = count ~ Intercept_count + elev_count + dist_count + space,
      family = "nzpoisson",
      data = counts_df[counts_df$present > 0, ],
      E = area)

presence_obs <- bru_obs(formula = present ~ Intercept_presence + elev_presence + dist_presence + space_copy,
  family = "binomial",
  data = counts_df)

fit_zap3 <- bru(
  cmp,
  presence_obs,
  pos_count_obs)


## -----------------------------------------------------------------------------
#| fig-cap: "Estimated expected counts for all four models"
pred_zip <- predict(
  fit_zip, 
  counts_df,
  ~ {
    pi <- zero_probability_parameter_for_zero_inflated_poisson_1
    lambda <- area * exp( distance + elevation + space + Intercept)
    expect <- (1-pi) * lambda
    variance <- (1-pi) * (lambda + pi * lambda^2)
    list(
      expect = expect
    )
  },n.samples = 2500)

pred_zap <- predict( fit_zap, counts_df,
  ~ {
    pi <- zero_probability_parameter_for_zero_inflated_poisson_0
    lambda <- area * exp( distance + elevation + space + Intercept)
    expect <- ((1-exp(-lambda))^(-1) * pi * lambda)
    list(
      expect = expect)
  },n.samples = 2500)

inv.logit = function(x) (exp(x)/(1+exp(x)))

pred_zap2 <- predict( fit_zap2, counts_df,
  ~ {
    pi <- inv.logit(Intercept_presence + elev_presence + dist_presence + space_presence)
    lambda <- area * exp( dist_count + elev_count + space_count + Intercept_count)
    expect <- ((1-exp(-lambda))^(-1) * pi * lambda)
    list(
      expect = expect)
  },n.samples = 2500)

pred_zap3 <- predict( fit_zap3, counts_df,
  ~ {
    pi <- inv.logit(Intercept_presence + elev_presence + dist_presence + space_copy)
    lambda <- area * exp( dist_count + elev_count + space + Intercept_count)
    expect <- ((1-exp(-lambda))^(-1) * pi * lambda)
    list(
      expect = expect)
  },n.samples = 2500)




  data.frame(x = st_coordinates(counts_df)[,1],
             y = st_coordinates(counts_df)[,2],
    zip = pred_zip$expect$mean,
         hurdle = pred_zap$expect$mean,
         hurdle2 = pred_zap2$expect$mean,
         hurdle3 = pred_zap3$expect$mean)  %>%
  pivot_longer(-c(x,y)) %>%
  ggplot() + geom_tile(aes(x,y, fill = value)) + facet_wrap(.~name) +
    theme_map + scale_fill_scico(direction = -1)








## -----------------------------------------------------------------------------
#| echo: true
#| eval: false

# bru_options_set(control.compute = list(dic = TRUE,
#                                        waic = TRUE,
#                                        mlik = TRUE,
#                                        cpo = TRUE))
# 


