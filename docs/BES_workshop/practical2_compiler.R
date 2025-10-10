## ----child="practicals\\spatial_data_types_areal.qmd"-------------------------

## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(CARBayesdata)

data(pollutionhealthdata)
data(GGHB.IZ)




## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)   
library(mapview)
library(sf)

# load some libraries to generate nice map plots
library(scico)


## -----------------------------------------------------------------------------
resp_cases <- merge(GGHB.IZ %>%
                      mutate(space = 1:dim(GGHB.IZ)[1]),
                             pollutionhealthdata, by = "IZ")%>%
  dplyr::filter(year == 2007) 


## -----------------------------------------------------------------------------
#| message: false
#| warning: false
resp_cases <- resp_cases %>% 
  mutate(SMR = observed/expected )


## -----------------------------------------------------------------------------
#| fig-width: 5
#| fig-height: 5
#| fig-align: center
#| code-fold: show
ggplot()+
  geom_sf(data=resp_cases,aes(fill=SMR))+
  scale_fill_scico(palette = "roma")


## -----------------------------------------------------------------------------
#| message: false
#| warning: false

library(spdep)

W.nb <- poly2nb(GGHB.IZ,queen = TRUE)
W.nb




## -----------------------------------------------------------------------------
library(spdep)
R <- nb2mat(W.nb, style = "B", zero.policy = TRUE)
diag = apply(R,1,sum)
Q = -R
diag(Q) = diag


## -----------------------------------------------------------------------------
#| echo: true
#| eval: false

# 
# cmp = ~ Intercept(1) + space(...) + iid(...)
# 
# formula = ...
# 
# 
# lik = bru_obs(formula = formula,
#               family = ...,
#               E = ...,
#               data = ...)
# 
# fit = bru(cmp, lik)
# 




## -----------------------------------------------------------------------------
#| echo: true
#| eval: false
# fit$summary.hyperpar




## -----------------------------------------------------------------------------
# produce predictions
pred = predict(fit, 
               resp_cases,
               ~data.frame(log_risk = Intercept + space,
                           risk = exp(Intercept + space),
                           cases = expected * exp(Intercept + space)),
               n.samples = 1000)

# plot the predictions

p1 = ggplot() + 
  geom_sf(data = pred$log_risk, aes(fill = mean)) + scale_fill_scico(direction = -1) +
  ggtitle("mean log risk")
p2 = ggplot() + 
  geom_sf(data = pred$log_risk, aes(fill = sd)) + scale_fill_scico(direction = -1) +
  ggtitle("sd log risk")

p3 = ggplot() +
  geom_sf(data = pred$risk, aes(fill = mean)) + scale_fill_scico(direction = -1) +
  ggtitle("mean  risk")

p4 = ggplot() + 
  geom_sf(data = pred$risk, aes(fill = sd)) +
  scale_fill_scico(direction = -1) +
  ggtitle("sd  risk")

p5 = ggplot() + geom_sf(data = pred$cases, aes(fill = mean)) + scale_fill_scico(direction = -1)+ 
  ggtitle("mean  expected counts")
p6 = ggplot() + geom_sf(data = pred$cases, aes(fill = sd)) + scale_fill_scico(direction = -1)+
  ggtitle("sd  expected counts")

p1 + p2 + p3 + p4 +p5 + p6 + plot_layout(ncol=2)



## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
pred$cases %>% ggplot() + geom_point(aes(observed, mean)) + 
  geom_errorbar(aes(observed, ymin = q0.025, ymax = q0.975)) +
  geom_abline(intercept = 0, slope = 1)



## -----------------------------------------------------------------------------
#| fig-align: center
#| fig-width: 4
#| fig-height: 4
#| code-fold: true
#| code-summary: "Click here to see the code"

# simulate 1000 realizations of E_i\lambda_i
expected_counts = generate(fit, resp_cases, 
                           ~ expected * exp(Intercept + space),
                           n.samples = 1000)


# simulate poisson data
aa = rpois(271*1000, lambda = as.vector(expected_counts))
sim_counts = matrix(aa, 271, 1000)

# summarise the samples with posterior means and quantiles
pred_counts = data.frame(observed = resp_cases$observed,
                         m = apply(sim_counts,1,mean),
                         q1 = apply(sim_counts,1,quantile, 0.025),
                         q2 = apply(sim_counts,1,quantile, 0.975),
                         vv = apply(sim_counts,1,var)
                         )
# Plot the observations against the predicted new counts and the predicted expected counts

ggplot() + 
  geom_point(data = pred_counts, aes(observed, m, color = "Pred_obs")) + 
  geom_errorbar(data = pred_counts, aes(observed, ymin = q1, ymax = q2, color = "Pred_obs")) +
  geom_point(data = pred$cases, aes(observed, mean, color = "Pred_means")) + 
  geom_errorbar(data = pred$cases, aes(observed, ymin = q0.025, ymax = q0.975, color = "Pred_means")) +
  
  geom_abline(intercept = 0, slope =1)




## ----child="practicals\\spatial_data_types_geostats.qmd"----------------------

## -----------------------------------------------------------------------------
#| message: false
#| warning: false
#| echo: false

library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)   
library(mapview)
library(sf)

# load some libraries to generate nice map plots
library(scico)



## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(sdmTMB)

pcod_df = sdmTMB::pcod %>% filter(year==2003)
qcs_grid = sdmTMB::qcs_grid



## -----------------------------------------------------------------------------
#| message: false
#| warning: false
#| 
pcod_sf =   st_as_sf(pcod_df, coords = c("lon","lat"), crs = 4326)
pcod_sf = st_transform(pcod_sf,
                       crs = "+proj=utm +zone=9 +datum=WGS84 +no_defs +type=crs +units=km" )


## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(terra)
depth_r <- rast(qcs_grid, type = "xyz")
crs(depth_r) <- crs(pcod_sf)


## -----------------------------------------------------------------------------
#| message: false
#| warning: false
#| code-fold: show
#| fig-align: center

library(tidyterra)
ggplot()+ 
  geom_spatraster(data=depth_r$depth)+
  geom_sf(data=pcod_sf,aes(color=factor(present))) +
    scale_color_manual(name="Occupancy status for the Pacific Cod",
                     values = c("black","orange"),
                     labels= c("Absence","Presence"))+
  scale_fill_scico(name = "Depth",
                   palette = "nuuk",
                   na.value = "transparent" ) + xlab("") + ylab("")



## -----------------------------------------------------------------------------
#| fig-align: center
#| fig-width: 5
#| fig-height: 5
mesh = fm_mesh_2d(loc = pcod_sf,           # Build the mesh
                  cutoff = 2,
                  max.edge = c(7,20),     # The largest allowed triangle edge length.
                  offset = c(5,50))       # The automatic extension distance

ggplot() + gg(mesh) +
  geom_sf(data= pcod_sf, aes(color = factor(present))) + 
  xlab("") + ylab("")


## -----------------------------------------------------------------------------
#| eval: false

# ?fm_mesh_2d


## -----------------------------------------------------------------------------
spde_model1 =  inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(.1, 0.5),
                                  prior.range = c(30, 0.5))
spde_model2 =  inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(10, 0.5),
                                  prior.range = c(1000, 0.5))
spde_model3 =  inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(100, 0.5))


## -----------------------------------------------------------------------------
cmp = ~ Intercept(1) + space(geometry, model = spde_model3)


## -----------------------------------------------------------------------------

formula = present ~ Intercept  + space

lik = bru_obs(formula = formula, 
              data = pcod_sf, 
              family = "binomial")



## -----------------------------------------------------------------------------
fit1 = bru(cmp,lik)




## -----------------------------------------------------------------------------
pxl = fm_pixels(mesh)


## -----------------------------------------------------------------------------
preds = predict(fit1, pxl, ~data.frame(spatial = space,
                                      total = Intercept + space))


## -----------------------------------------------------------------------------
#| fig-align: center
#| code-fold: show

ggplot() + geom_sf(data = preds$spatial,aes(color = mean)) + 
  scale_color_scico() +
  ggtitle("Posterior mean") +
ggplot() + geom_sf(data = preds$spatial,aes(color = sd)) + 
  scale_color_scico() +
  ggtitle("Posterior sd")




## -----------------------------------------------------------------------------
#| fig-align: center
#| fig-width: 4
#| fig-height: 4
# create the grouped variable
depth_r$depth_group = inla.group(values(depth_r$depth_scaled))

# run the model
cmp = ~ Intercept(1) + space(geometry, model = spde_model3) +
        covariate(depth_r$depth_group, model = "rw2")

formula = present ~ Intercept  + space + covariate

lik = bru_obs(formula = formula, 
              data = pcod_sf, 
              family = "binomial")


fit3 = bru(cmp, lik)

# plot the estimated effect of depth

fit3$summary.random$covariate %>% 
  ggplot() + geom_line(aes(ID,mean)) + 
             geom_ribbon(aes(ID,
                             ymin = `0.025quant`, 
                             ymax = `0.975quant`),
                         alpha = 0.5)


## -----------------------------------------------------------------------------
pxl1 = data.frame(crds(depth_r), 
                  as.data.frame(depth_r$depth)) %>% 
       filter(!is.na(depth)) %>%
st_as_sf(coords = c("x","y")) %>% dplyr::select(-depth)


## -----------------------------------------------------------------------------
inv_logit = function(x) (1+exp(-x))^(-1)




## ----child="practicals\\spatial_data_types_PP.qmd"----------------------------

## ----test---------------------------------------------------------------------
#| message: false
#| warning: false
#| echo: false
library(webexercises)
# For plotting
library(mapview)
library(ggplot2)
library(scico) # for colouring palettes
library(patchwork)
library(terra)
library(sf)
# Data manipulation
library(dplyr)
library(fmesher)
library(INLA)
library(inlabru)
library(viridis)
library(gt)


## -----------------------------------------------------------------------------
# load data files (set the correct path to where the data are stored)
load("datasets/ringlett_CNP.RData")
elev_rast = terra::rast("datasets/elev_CNP.tiff")


## -----------------------------------------------------------------------------
shp_SGC <-  st_read("datasets/SG_CairngormsNationalPark/SG_CairngormsNationalPark_2010.shp",quiet =T)


## -----------------------------------------------------------------------------
shp_SGC <- shp_SGC %>% st_transform(crs = 27700)
shp_SGC <- st_transform(shp_SGC,gsub("units=m","units=km",st_crs(shp_SGC)$proj4string)) 



## -----------------------------------------------------------------------------
#| code-fold: true
#| fig-align: center
#| fig-width: 5
#| fig-height: 5
ggplot()+
  geom_sf(data=shp_SGC)+
  geom_sf(data=ringlett_CNP,size=2)




## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
# mesh options
max.edge = 2
bound.outer = 10

mesh <- fm_mesh_2d(boundary = shp_SGC,
                   max.edge = c(1,3)*max.edge,
                   offset=c(max.edge, bound.outer),
                   cutoff = 1.5,
                   crs=st_crs(shp_SGC))

ggplot() + gg(mesh) + geom_sf(data=ringlett_CNP)


## -----------------------------------------------------------------------------
st_area(shp_SGC)


## -----------------------------------------------------------------------------
matern <- inla.spde2.pcmatern(mesh,
                              prior.range = c(20, 0.5), # P(range < 20) = 0.5
                              prior.sigma = c(1, 0.5))  # P(sigma > 1) = 0.5


## -----------------------------------------------------------------------------
#| message: false
#| warning: false
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
library(stars)

# Standardize the raster 
elev_rast <- scale(elev_rast)
# Extend raster ext by a factor of 1.5
re <- extend(elev_rast, ext(elev_rast)*1.5)
# Convert to spatial object
re_df <- re %>% stars::st_as_stars() %>%  st_as_sf(na.rm=F)
# fill in missing values using the original raster 
re_df$GBR_elv_msk <- bru_fill_missing(elev_rast,re_df,re_df$GBR_elv_msk)
# rasterize
elev_rast_p <- stars::st_rasterize(re_df) %>% rast()
# visualize
plot(elev_rast_p)



## -----------------------------------------------------------------------------
cmp_lgcp <-  geometry ~  Intercept(1)  + 
  elevation(elev_rast_p, model = "linear") +
  grf(geometry, model = matern)


## -----------------------------------------------------------------------------
formula = geometry ~ Intercept  + elevation + grf


## -----------------------------------------------------------------------------
lik = bru_obs(formula = formula, 
              data = ringlett_CNP, 
              family = "cp",
              samplers = shp_SGC,
              domain = list(geometry = mesh))



## -----------------------------------------------------------------------------
fit_lgcp = bru(cmp_lgcp,lik)


## -----------------------------------------------------------------------------
#| eval: false

# summary(fit_lgcp)




## -----------------------------------------------------------------------------
pxl1 = data.frame(crds(elev_rast), 
                  as.data.frame(elev_rast$GBR_elv_msk)) %>% 
       filter(!is.na(GBR_elv_msk)) %>%
st_as_sf(coords = c("x","y")) %>%
  dplyr::select(-GBR_elv_msk)




## -----------------------------------------------------------------------------
#| eval: false

# lgcp_pred <- predict(
#   fit_lgcp,
#   pxl1,
#   ~ data.frame(
#     lambda = exp(Intercept + elevation + grf), # intensity
#     loglambda = Intercept + elevation +grf,  #log-intensity
#     GF = grf # matern field
#   )
# )
# 
# # predicted log intensity
# ggplot() + gg(lgcp_pred$loglambda, geom = "tile")
# # standard deviation of the predicted log intensity
# ggplot() + gg(lgcp_pred$loglambda, geom = "tile",aes(fill=sd))
# # predicted intensity
# ggplot() +  gg(lgcp_pred$lambda, geom = "tile")
# # spatial field
# ggplot() +  gg(lgcp_pred$GF, geom = "tile")


