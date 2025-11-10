## ----child="practicals/zero_inflated_ex.qmd"----------------------------------

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
  scale_fill_scico(direction = -1) + geom_sf(data = boundary, alpha = 0) + theme_map+ theme(legend.position = "none")
p2 = ggplot() + geom_spatraster(data = dist_cov) + ggtitle("Distance to water") +
  scale_fill_scico(direction = -1) + geom_sf(data = boundary, alpha = 0) +
  theme_map + theme(legend.position = "none")
  
p1 + p2


## -----------------------------------------------------------------------------
#| fig-cap: "Counts of gorilla nests"
#| out-width: "70%"
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
#| out-width: 80%
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




p =   data.frame(x = st_coordinates(counts_df)[,1],
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


