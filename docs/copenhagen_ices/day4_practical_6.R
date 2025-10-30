## ----child="practicals/LM_resids.qmd"-----------------------------------------

## -----------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false

library(webexercises)



## -----------------------------------------------------------------------------
#| warning: false
#| message: false

library(dplyr)
library(tidyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     


## -----------------------------------------------------------------------------
#| code-fold: show
beta = c(2,0.5)
sd_error = 0.1

n = 100
x = rnorm(n)
y = beta[1] + beta[2] * x + rnorm(n, sd = sd_error)

df = data.frame(y = y, x = x)  



## -----------------------------------------------------------------------------
# Model components
cmp =  ~ -1 + beta_0(1) + beta_1(x, model = "linear")
# Linear predictor
formula = y ~ Intercept + beta_1
# Observational model likelihood
lik =  bru_obs(formula = y ~.,
            family = "gaussian",
            data = df)
# Fit the Model
fit.lm = bru(cmp, lik)


## -----------------------------------------------------------------------------
res_samples <- predict(
  fit.lm,         # the fitted model
  df,             # the original data set
  ~ data.frame(   
    res = y-(beta_0 + beta_1)  # compute the residuals
  ),
  n.samples = 1000   # draw 1000 samples
)



## -----------------------------------------------------------------------------
#| fig-width: 6
#| fig-height: 4
#| fig-align: center
#| fig-cap: "Bayesian residual plots: the left panel is the residual index plot; the right panel is the plot of the residual versus the covariate x"
#| code-fold: true
#| code-summary: "Residuals checks for Linear Model"

ggplot(res_samples,aes(y=mean,x=1:100))+geom_point() +
ggplot(res_samples,aes(y=mean,x=x))+geom_point()



## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
#| code-fold: true
#| code-summary: "QQPlot for Linear Model"

arrange(res_samples, mean) %>%
  mutate(theortical_quantiles = qnorm(1:100 / (1+100))) %>%
  ggplot(aes(x=theortical_quantiles,y= mean)) + 
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "grey70")+
  geom_abline(intercept = mean(res_samples$mean),
              slope = sd(res_samples$mean)) +
  geom_point() +
  labs(x = "Theoretical Quantiles (Normal)",
       y= "Sample Quantiles (Residuals)") 


## -----------------------------------------------------------------------------
#| warning: false
#| message: false

samples =  generate(fit.lm, df,
  formula = ~ {
    mu <- (beta_0 + beta_1)
    sd <- sqrt(1 / Precision_for_the_Gaussian_observations)
    rnorm(100, mean = mu, sd = sd)
  },
  n.samples = 500
) 


## -----------------------------------------------------------------------------
# Tidy format for plotting
samples_long = data.frame(samples) %>% 
  mutate(id = 1:100) %>% # i-th observation
  pivot_longer(-id)

# compute the mean and quantiles for the samples
draws_summaries = data.frame(mean_samples = apply(samples,1,mean),
q25 = apply(samples,1,function(x)quantile(x,0.025)),  
q975 = apply(samples,1,function(x)quantile(x,0.975)),
observations = df$y)  

p1 = ggplot() + geom_density(data = samples_long, 
                        aes(value, group = name),  color = "#E69F00") +
  geom_density(data = df, aes(y))  +
  xlab("") + ylab("") 

p2 = ggplot(draws_summaries,aes(y=mean_samples,x=observations))+
  geom_errorbar(aes(ymin = q25,
                   ymax = q975), 
               alpha = 0.5, color = "grey50")+
geom_point()+geom_abline(slope = 1,intercept = 0,lty=2)+labs()

p1 +p2



## ----child="practicals/GLM_scores_ex.qmd"-------------------------------------

## -----------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false

library(webexercises)



## -----------------------------------------------------------------------------
#| warning: false
#| message: false

library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     



## -----------------------------------------------------------------------------
crabs <- read.csv("datasets/crabs.csv")

# conditional means and variances
crabs %>%
  summarise( Mean = mean(satell ),
             Variance = var(satell),
                     .by = color)


## -----------------------------------------------------------------------------

crabs_df = model.matrix( ~  color , crabs) %>%
  as.data.frame() %>%
  select(-1) %>%        # drop intercept
  bind_cols(crabs) %>%  # append to original data
  select(-color)        # remove original color categorical variable



## -----------------------------------------------------------------------------

cmp =  ~ -1 + beta0(1) +  colordarker +
       colorlight + colormedium +
       w(weight, model = "linear")

lik =  bru_obs(formula = satell ~.,
            family = "poisson",
            data = crabs_df)

fit_pois = bru(cmp, lik)

summary(fit_pois)



## -----------------------------------------------------------------------------
bru_options_set(control.compute = list(cpo = TRUE))

fit_pois = bru(cmp, lik)




## -----------------------------------------------------------------------------
#| eval: false

# 
# fit_pois$cpo$pit %>%
#   hist(main = "Histogram of PIT values")
# 
# qqplot(qunif(ppoints(length(fit_pois$cpo$pit))),
#        fit_pois$cpo$pit,
#        main = "Q-Q plot for Unif(0,1)",
#        xlab = "Theoretical Quantiles",
#        ylab = "Sample Quantiles")
# 
# qqline(fit_pois$cpo$pit,
#        distribution = function(p) qunif(p),
#        prob = c(0.1, 0.9))



## ----child="practicals/LGOCV_ex.qmd"------------------------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(inlabru) 
library(sf)
library(terra)


# load some libraries to generate nice map plots
library(scico)
library(ggplot2)
library(patchwork)
library(mapview)
library(tidyterra)



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

depth_r <- rast(qcs_grid, type = "xyz")
crs(depth_r) <- crs(pcod_sf)





## -----------------------------------------------------------------------------
mesh = fm_mesh_2d(loc = pcod_sf,         
                  max.edge = c(10,20),     
                  offset = c(5,50))   

spde_model =  inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(100, 0.5))




## -----------------------------------------------------------------------------
# Model 1 components
cmp_1 = ~ Intercept(1) + space(geometry, model = spde_model)

# Model 1 linear predictor
formula_1 = present ~ Intercept  + space

# Model 1 observational model
lik_1 = bru_obs(formula = formula_1, 
              data = pcod_sf, 
              family = "binomial")
# fit Model 1
fit_spatial = bru(cmp_1,lik_1)


## -----------------------------------------------------------------------------
# Model 2 components
cmp_2 = ~ Intercept(1) + 
  space(geometry, model = spde_model)+
  covariate(depth_r$depth_scaled, model = "linear") 

# Model 2 linear predictor
formula_2 = present ~ Intercept  + covariate + space

# Model 2 observational model
lik_2 = bru_obs(formula = formula_2, 
              data = pcod_sf, 
              family = "binomial")

# Fit Model 2
fit_depth_linear = bru(cmp_2,lik_2)


## -----------------------------------------------------------------------------
# create the grouped variable
depth_r$depth_group = inla.group(values(depth_r$depth_scaled))

# Model 3 components
cmp_3 = ~ Intercept(1) + 
  space(geometry, model = spde_model)+
  covariate(depth_r$depth_group, model = "rw2")

# Model 3 linear predictor
formula_3 = present ~ Intercept  + covariate + space

# Model 2 observational model
lik_3 = bru_obs(formula = formula_2, 
              data = pcod_sf, 
              family = "binomial")

# Fit Model 2
fit_depth_smooth = bru(cmp_3,lik_3)


## -----------------------------------------------------------------------------
# create buffer of size 25 centred at each site
buffer <- st_buffer(pcod_sf, dist = 25)

# Lists of the indexes of the leave-out-group for each observation i
Ii <- st_intersects(pcod_sf,buffer)




## -----------------------------------------------------------------------------
lgocv_depth_smooth = inla.group.cv(result = fit_depth_smooth,
                                   groups = Ii)
log_depth_smooth = mean(log(lgocv_depth_smooth$cv),
                        na.rm=T)


## -----------------------------------------------------------------------------
#| code-fold: false
#| warning: false
lgocv_spatial = inla.group.cv(result = fit_spatial, 
                              group.cv = lgocv_depth_smooth)

lgocv_depth_linear = inla.group.cv(result = fit_depth_linear,
                                   group.cv = lgocv_depth_smooth)

log_score_spatial<- mean(log(lgocv_spatial$cv),na.rm=T)
log_depth_linear <-mean(log(lgocv_depth_linear$cv),na.rm=T)
log_depth_smooth <-mean(log(lgocv_depth_smooth$cv),na.rm=T)



## -----------------------------------------------------------------------------
#| eval: false
# data.frame(logspat= log_score_spatial,
#            logdepthl = log_depth_linear,
#            logdepths = log_depth_smooth )




## -----------------------------------------------------------------------------
pcod_spat = sdmTMB::pcod %>%
  filter(year %in% 2003:2011) %>%
  st_as_sf( coords = c("lon","lat"), crs = 4326) %>%
  st_transform(., crs = "+proj=utm +zone=9 +datum=WGS84 +no_defs +type=crs +units=km" ) %>%
   mutate(time_idx = match(year, c(2003, 2004, 2005, 2007, 2009, 2011, 2013, 2015, 2017)),
         id = 1:nrow(.)) # Observation id for CV


## -----------------------------------------------------------------------------
bru_options_set(control.compute = list(waic = TRUE,dic= TRUE,mlik = TRUE))


## -----------------------------------------------------------------------------
#| message: false
#| warning: false 
# Model components
cmp_spat = ~ Intercept(1) + 
  covariate(depth_r$depth_group, model = "rw2")+
  trend(time_idx, model = "iid")+
  space(geometry, model = spde_model)

# Linear predictor
formula_spat = present ~ Intercept  + covariate  + trend + space

# Observational model
lik_spat = bru_obs(formula = formula_spat, 
              data = pcod_spat, 
              family = "binomial")

# Fit Model 
fit_spat = bru(cmp_spat,lik_spat)



## -----------------------------------------------------------------------------
#| eval: false
# # PC prior for AR(1) correlation parameter
# h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
# 
# # Model components
# cmp_spat_ar1 = ~ Intercept(1) +
#   covariate(depth_r$depth_group, model = "rw2")+
#   space_time(geometry,
#         group = time_idx ,
#         model = spde_model,
#         control.group = list(model = 'ar1',hyper = h.spec))
# 
# # Linear predictor
# formula_spat_ar1 = present ~ Intercept  + covariate  + space_time
# 
# # Observational model
# lik_spat_ar1 = bru_obs(formula = formula_spat_ar1,
#               data = pcod_spat,
#               family = "binomial")
# 
# # Fit Model
# fit_spat_ar1 = bru(cmp_spat_ar1,lik_spat_ar1)
# 


## -----------------------------------------------------------------------------

# create buffer of size 25  centred at each site

buffer_25 <- st_buffer(pcod_spat, dist = 25) 


# empty lists to include the indexes of the leave-out-group for each observation i
I_i <- list()

# loop though each observation and store the leave-out-group based on the buffer
for( i in 1:nrow(pcod_spat)){
  
  # Temporal filtering of data within a 2 years of span of  observation i
  df_sf_subset <- pcod_spat %>% 
    filter( between(time_idx,left = pcod_spat$time_idx[i]-2, 
                    right = pcod_spat$time_idx[i]+2)) 
  # Spatial filtering of the observations that are within the buffer of the ith observation
  Buffer_i <-df_sf_subset %>% st_intersects(buffer_25[i,],sparse = FALSE) %>% # identify 
    unlist()
  
  # obtain the indexes of the leave out group
  I_i[[i]] <-  df_sf_subset[Buffer_i,] %>%  pull(id)
  
}





## -----------------------------------------------------------------------------
#| eval: false
# lgocv_spat_ar1 <- inla.group.cv(result = fit_spat_ar1, groups  = I_i)
# logscore_spat_ar1 = mean(log(lgocv_spat_ar1$cv),na.rm=T)
# 
# 
# lgocv_spat <- inla.group.cv(result = fit_spat, group.cv  = lgocv_spat_ar1)
# logscore_spat = mean(log(lgocv_spat$cv),na.rm=T)
# 
# 


## -----------------------------------------------------------------------------
#| eval: false
# table = data.frame( DIC = c(fit_spat_ar1$dic$dic, fit_spat$dic$dic),
#                     WAIC = c(fit_spat_ar1$waic$waic, fit_spat$waic$waic),
#                     mlik = c(fit_spat_ar1$mlik[1,1],fit_spat$mlik[1,1]),
#                     LGOCV = c(logscore_spat_ar1,logscore_spat))
# 
#  rownames(table) = c("Model 2 - spatiotemporal field",
#                      "Model 1 - time iid effect")
# 


