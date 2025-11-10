## ----child="practicals/Multiple_likelihood.qmd"-------------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(inlabru) 
library(sf)
library(terra)
library(tidyverse)

# load some libraries to generate nice map plots
library(scico)
library(ggplot2)
library(patchwork)



## -----------------------------------------------------------------------------

N = 200
x =  runif(N)
df = data.frame(idx = 1:N,
                x = x)

# simulate data
df = df %>% 
  mutate(y_gaus = rnorm(N, mean = 1 + 1.5 * x), sd = 0.5) %>%
  mutate(y_pois = rpois(N, lambda  = exp( -1 + 1.5 * x))) 

# plot the data
df %>% ggplot() + 
  geom_point(aes(x, y_gaus, color = "Gaussian")) +
  geom_point(aes(x, y_pois, color = "Poisson")) 



## -----------------------------------------------------------------------------
cmp = ~ -1 + 
  Intercept_gaus(1) + 
  Intercept_pois(1) +
  covariate(x, model = "linear") 


## -----------------------------------------------------------------------------
lik_gaus = bru_obs(formula = y_gaus ~ Intercept_gaus + covariate,
                    data = df)

lik_pois = bru_obs(formula = y_pois ~ Intercept_pois + covariate,
                    data = df,
                   family = "poisson")




## -----------------------------------------------------------------------------
rMatern <- function(n, coords, sigma=1, range, 
                    kappa = sqrt(8*nu)/range, 
                    variance = sigma^2, 
                    nu=1) {
  m <- as.matrix(dist(coords))
  m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  diag(m) <- 1
  return(drop(crossprod(chol(variance*m),
                        matrix(rnorm(nrow(coords)*n), ncol=n))))
}


## -----------------------------------------------------------------------------
# Intercept on reparametrized model
beta <- c(-5, 3) 
# Random field marginal variances for omega1 and omega2:
m.var <- c(0.5, 0.4) 
# GRF range parameters for omega1 and omega2:
range <- c(4, 6)
# Copy parameters: reparameterization of coregionalization 
# parameters
lambda <- c(0.7) 
# Standard deviations of error terms
e.sd <- c(0.3, 0.2)



## -----------------------------------------------------------------------------
# define the area of interest
poly_geom = st_polygon(list(cbind(c(0,10,10,0,0), c(0,0,5,5,0)) ))
# Wrap it in an sfc (simple feature collection)
poly_sfc <- st_sfc(poly_geom)
# Now create the sf object
border <- st_sf(id = 1, geometry = poly_sfc)



# how many observation we have
n1 <- 200
n2 <- 150
n_common = 50

# simulate observation locations

loc_common = st_sf(geometry = st_sample(border, n_common))
loc_only1 = st_sf(geometry = st_sample(border, n1-n_common))
loc_only2 = st_sf(geometry = st_sample(border, n2-n_common))



# simulate the two gaussian field at the locations
z1 <- rMatern(1, st_coordinates(rbind( loc_common,loc_only1, loc_only2)), range = range[1],
                  sigma = sqrt(m.var[1]))

z2 <- rMatern(1, st_coordinates(rbind(loc_common, loc_only2)), range = range[2],
                  sigma = sqrt(m.var[2]))


## Create data.frame
loc1 = rbind( loc_common, loc_only1)
loc2 = rbind( loc_common, loc_only2)

df1 =  loc1 %>% mutate(z1 = z1[1:n1])
df2 =  loc2 %>% mutate(z1 = z1[-c(1:(n1-n_common))], z2 =z2)


## create the linear predictors

df1  = df1 %>%
  mutate(eta1 = beta[1] + z1)

df2  = df2 %>%
  mutate(eta2 = beta[2] + lambda * z1 + z2)


# simulate data by addint the obervation noise

df1  = df1 %>%
  mutate(y = rnorm(n1, mean = eta1, sd = e.sd[1]))

df2  = df2 %>%
  mutate(y = rnorm(n2, mean = eta2, sd = e.sd[1]))


## ----out.width="95%"----------------------------------------------------------
p1 = ggplot(data = df1) + geom_sf(aes(color = z1)) 
p2 = ggplot(data = df2) + geom_sf(aes(color = z2)) 
p1+p2+plot_layout(ncol = 1)


## -----------------------------------------------------------------------------
mesh <-  fm_mesh_2d(loc = rbind(loc1, loc2), 
                   boundary = border,
                     max.edge = c(0.5, 1.5), 
                     offset = c(0.1, 2.5), 
                     cutoff = 0.1)



## ----echo = F-----------------------------------------------------------------
ggplot() + 
  gg(mesh) + 
  geom_sf(data =  df1,  size = 2, aes(color = "data 1")) +
    geom_sf(data =  df2, aes(color = "data 2")) + xlab("") + ylab("")




## -----------------------------------------------------------------------------
cmp = ~ -1 +  Intercept1(1) + Intercept2(1) +
  omega1(geometry, model = spde) +
  omega1_copy(geometry, copy = "omega1", fixed = FALSE) +
  omega2(geometry, model = spde)



