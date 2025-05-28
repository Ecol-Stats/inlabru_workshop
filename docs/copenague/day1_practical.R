## ----child="C:\\Users\\jb538u\\OneDrive - University of Glasgow\\Documents\\inlabru_workshop2025\\inlabru_workshop\\docs\\practicals\\test.qmd"----

## --------------------------------------------------------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     
# load some libraries to generate nice map plots
library(scico)


## --------------------------------------------------------------------------------------------------------------------------
#| code-fold: show

beta = c(1,1)
sd_error = 1

n = 100
x = rnorm(n)
y = beta[1] + beta[2] * x + rnorm(n, sd = sd_error)

df = data.frame(y = y, x = x)  



## --------------------------------------------------------------------------------------------------------------------------
cmp =  ~ Intercept(1) + beta_1(x, model = "linear")




## --------------------------------------------------------------------------------------------------------------------------
#| eval: false
# formula = y ~ Intercept + beta_1


## --------------------------------------------------------------------------------------------------------------------------
lik =  bru_obs(formula = y ~.,
            family = "gaussian",
            data = df)


## --------------------------------------------------------------------------------------------------------------------------
fit.lm = bru(cmp, lik)


## --------------------------------------------------------------------------------------------------------------------------
summary(fit.lm)


## --------------------------------------------------------------------------------------------------------------------------
plot(fit.lm, "Intercept")






## --------------------------------------------------------------------------------------------------------------------------
new_data = data.frame(x = c(df$x, runif(10)),
                      y = c(df$y, rep(NA,10)))
pred = predict(fit.lm, new_data, ~ Intercept + beta_1)


## --------------------------------------------------------------------------------------------------------------------------
#| code-fold: true
#| fig-cap: Data and 95% credible intervals
#| echo: false
#| message: false
#| warning: false

pred %>% ggplot() + 
  geom_point(aes(x,y), alpha = 0.3) +
  geom_line(aes(x,mean)) +
  geom_line(aes(x, q0.025), linetype = "dashed")+
  geom_line(aes(x, q0.975), linetype = "dashed")+
  xlab("Covariate") + ylab("Observations")


## --------------------------------------------------------------------------------------------------------------------------
#| code-fold: show
#| eval: false
#| 
# pred %>% ggplot() +
#   geom_point(aes(x,y), alpha = 0.3) +
#   geom_line(aes(x,mean)) +
#   geom_line(aes(x, q0.025), linetype = "dashed")+
#   geom_line(aes(x, q0.975), linetype = "dashed")+
#   xlab("Covariate") + ylab("Observations")



## ----child="C:\\Users\\jb538u\\OneDrive - University of Glasgow\\Documents\\inlabru_workshop2025\\inlabru_workshop\\docs\\practicals\\LMM_ex.qmd"----

## ----load_libs_lmm, message=FALSE, warning=FALSE---------------------------------------------------------------------------
library(inlabru)
library(INLA)
library(ggplot2)
library(dplyr)


## ----generate_data_lmm-----------------------------------------------------------------------------------------------------

beta = c(1.5,1)
sd_error = 1
tau_group = 1

n = 100
n.groups = 5
x = rnorm(n)
v = rnorm(n.groups, sd = tau_group^{-1/2})
y = beta[1] + beta[2] * x + rnorm(n, sd = sd_error) +
  rep(v, each = 20)

df = data.frame(y = y, x = x, j = rep(1:5, each = 20))  


## ----plot_data_lmm---------------------------------------------------------------------------------------------------------
#| code-fold: true
#| fig-cap: Data for the linear mixed model example with 5 groups
df$jfac = as.factor(df$j)
ggplot(df) +
  geom_point(aes(x = x, colour = jfac, y = y)) +
  theme_classic() +
  scale_colour_discrete("Group")



## ----define_components_lmm-------------------------------------------------------------------------------------------------
# Define model components
cmp =  ~ Intercept(1) + beta_1(x, model = "linear") +
  v(j, model = "iid")


## ----define_likelihood_lmm-------------------------------------------------------------------------------------------------
# Construct likelihood
lik =  like(formula = y ~.,
            family = "gaussian",
            data = df)


## ----fit_model_lmm---------------------------------------------------------------------------------------------------------
fit = bru(cmp, lik)
summary(fit)


