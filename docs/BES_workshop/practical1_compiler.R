## ----child="practicals\\LM_ex.qmd"--------------------------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false
#| code-summary: "Load libraries"

library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     
# load some libraries to generate nice plots
library(scico)


## -----------------------------------------------------------------------------
#| code-fold: show
#| code-summary: "Simulate Data from a LM"

beta = c(2,0.5)
sd_error = 0.1

n = 100
x = rnorm(n)
y = beta[1] + beta[2] * x + rnorm(n, sd = sd_error)

df = data.frame(y = y, x = x)  



## -----------------------------------------------------------------------------
#| code-summary: "Define LM components"
cmp =  ~ -1 + beta_0(1) + beta_1(x, model = "linear")




## -----------------------------------------------------------------------------
#| eval: false
#| code-summary: "Define LM formula"
# formula = y ~ beta_0 + beta_1


## -----------------------------------------------------------------------------
#| code-summary: "Define Observational model"
lik =  bru_obs(formula = y ~.,
            family = "gaussian",
            data = df)


## -----------------------------------------------------------------------------
#| code-summary: "Fit LM in `inlabru`"
fit.lm = bru(cmp, lik)


## -----------------------------------------------------------------------------
#| code-summary: "Model summaries"
#| collapse: true
summary(fit.lm)


## -----------------------------------------------------------------------------
new_data = data.frame(x = c(df$x, runif(10)),
                      y = c(df$y, rep(NA,10)))
pred = predict(fit.lm, new_data, ~ beta_0 + beta_1,
               n.samples = 1000)





## -----------------------------------------------------------------------------
#| code-fold: true
#| fig-cap: Data and 95% credible intervals
#| echo: false
#| message: false
#| warning: false
#| fig-align: center
#| fig-width: 4
#| fig-height: 4

pred %>% ggplot() + 
  geom_point(aes(x,y), alpha = 0.3) +
  geom_line(aes(x,mean)) +
  geom_line(aes(x, q0.025), linetype = "dashed")+
  geom_line(aes(x, q0.975), linetype = "dashed")+
  xlab("Covariate") + ylab("Observations")






## -----------------------------------------------------------------------------
#| eval: true
#| purl: true 

inla.priors.used(fit.lm)






## -----------------------------------------------------------------------------

# First we define the logGamma (0.01,0.01) prior 

prec.tau <- list(prec = list(prior = "loggamma",   # prior name
                             param = c(0.01, 0.01))) # prior values

lik2 =  bru_obs(formula = y ~.,
                family = "gaussian",
                data = df,
                control.family = list(hyper = prec.tau))

fit.lm2 = bru(cmp2, lik2) 



## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-align: center
#| fig-height: 4
#| code-fold: show

plot(fit.lm, "beta_0")



## ----child="practicals\\LMM_ex.qmd"-------------------------------------------

## -----------------------------------------------------------------------------
#| code-summary: "Simulate data from a LMM"
#| 
set.seed(12)
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


## ----plot_data_lmm------------------------------------------------------------
#| code-fold: true
#| fig-cap: Data for the linear mixed model example with 5 groups
#| fig-align: center
#| fig-width: 4
#| fig-height: 4

ggplot(df) +
  geom_point(aes(x = x, colour = factor(j), y = y)) +
  theme_classic() +
  scale_colour_discrete("Group")



## ----define_components_lmm----------------------------------------------------
# Define model components
cmp =  ~ -1 + beta_0(1) + beta_1(x, model = "linear") +
  u(j, model = "iid")


## ----define_likelihood_lmm----------------------------------------------------
# Construct likelihood
lik =  like(formula = y ~.,
            family = "gaussian",
            data = df)


## -----------------------------------------------------------------------------
#| collapse: true
#| code-summary: "Fit a LMM in inlabru"
fit = bru(cmp, lik)
summary(fit)


## -----------------------------------------------------------------------------
#| code-fold: true
#| code-summary: "LMM fitted values"
#| fig-align: center
#| fig-width: 4
#| fig-height: 4

# New data
xpred = seq(range(x)[1], range(x)[2], length.out = 100)
j = 1:n.groups
pred_data = expand.grid(x = xpred, j = j)
pred = predict(fit, pred_data, formula = ~ beta_0 + beta_1 + u) 


pred %>%
  ggplot(aes(x=x,y=mean,color=factor(j)))+
  geom_line()+
  geom_ribbon(aes(x,ymin = q0.025, ymax= q0.975,fill=factor(j)), alpha = 0.5) + 
  geom_point(data=df,aes(x=x,y=y,colour=factor(j)))+
  facet_wrap(~j)




## ----child="practicals\\GLM_ex.qmd"-------------------------------------------

## -----------------------------------------------------------------------------
#| code-summary: "Simulate Data from a GLM"
set.seed(123)
n = 100
beta = c(1,1)
x = rnorm(n)
lambda = exp(beta[1] + beta[2] * x)
y = rpois(n, lambda  = lambda)
df = data.frame(y = y, x = x)  


## -----------------------------------------------------------------------------
#| code-summary: "GLM components"

cmp =  ~ -1 + beta_0(1) + beta_1(x, model = "linear")


## -----------------------------------------------------------------------------
#| code-summary: "GLM likelihood"

lik =  bru_obs(formula = y ~.,
            family = "poisson",
            data = df)


## -----------------------------------------------------------------------------
#| code-summary: "Fit a GLM"
fit_glm = bru(cmp, lik)


## -----------------------------------------------------------------------------
#| code-summary: "GLM summaries"
summary(fit_glm)


## ----get_predictions_glm------------------------------------------------------
#| code-summary: "Predcited values for Poisson GLM"

# Define new data, set to NA the values for prediction

new_data = data.frame(x = c(df$x, runif(10)),
                      y = c(df$y, rep(NA,10)))

# Define predictor formula
pred_fml <- ~ exp(beta_0 + beta_1)

# Generate predictions
pred_glm <- predict(fit_glm, new_data, pred_fml)


## -----------------------------------------------------------------------------
#| echo: false
#| code-summary: "Plot GLM predicted values"
#| fig-cap: Data and 95% credible intervals
#| fig-align: center
#| fig-width: 4
#| fig-height: 4
#| warning: false

pred_glm %>% ggplot() + 
  geom_point(aes(x,y), alpha = 0.3) +
  geom_line(aes(x,mean)) +
    geom_ribbon(aes(x = x, ymax = q0.975, ymin = q0.025),fill = "tomato", alpha = 0.3)+
  geom_line(aes(x, q0.025), linetype = "dashed")+
  geom_line(aes(x, q0.975), linetype = "dashed")+
  xlab("Covariate") + ylab("Observations (counts)")




## -----------------------------------------------------------------------------
#| code-summary: "GLM Task"
set.seed(123)
n = 100
alpha = c(0.5,1.5)
w = rnorm(n)
psi = plogis(alpha[1] + alpha[2] * w)
y = rbinom(n = n, size = 1, prob =  psi) # set size = 1 to draw binary observations
df_logis = data.frame(y = y, w = w)  



## ----child="practicals\\GAM_ex.qmd"-------------------------------------------

## -----------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false
library(webexercises)
library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     
# load some libraries to generate nice map plots
library(scico)


## -----------------------------------------------------------------------------
#| code-summary: "Simulate GAM Data"
n = 100
x = rnorm(n)
eta = (1 + cos(x))
y = rnorm(n, mean =  eta, sd = 0.5)

df = data.frame(y = y, 
                x_smooth = inla.group(x)) # equidistant x's 



## ----define_components_gam----------------------------------------------------
cmp =  ~ Intercept(1) + 
  smooth(x_smooth, model = "rw1")


## ----define_likelihood_gam----------------------------------------------------
lik =  bru_obs(formula = y ~.,
            family = "gaussian",
            data = df)


## ----run_model_gam------------------------------------------------------------
fit = bru(cmp, lik)
fit$summary.fixed




## -----------------------------------------------------------------------------
#| eval: false
#| 
# data.frame(fit$summary.random$smooth) %>%
#   ggplot() +
#   geom_ribbon(aes(ID,ymin = X0.025quant, ymax= X0.975quant), alpha = 0.5) +
#   geom_line(aes(ID,mean)) +
#   xlab("covariate") + ylab("")


## ----get_predictions_gam------------------------------------------------------
pred = predict(fit, df, ~ (Intercept + smooth))


## ----plot_gam-----------------------------------------------------------------
#| code-fold: true
#| fig-cap: Data and 95% credible intervals
pred %>% ggplot() + 
  geom_point(aes(x_smooth,y), alpha = 0.3) +
  geom_line(aes(x_smooth,1+cos(x_smooth)),col=2)+
  geom_line(aes(x_smooth,mean)) +
  geom_line(aes(x_smooth, q0.025), linetype = "dashed")+
  geom_line(aes(x_smooth, q0.975), linetype = "dashed")+
  xlab("Covariate") + ylab("Observations")


