## ----child= "C:\\Users\\jb538u\\OneDrive - University of Glasgow\\Documents\\inlabru_workshop2025\\inlabru_workshop\\docs\\practicals\\LM_ex.qmd"----

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



## ----child= "C:\\Users\\jb538u\\OneDrive - University of Glasgow\\Documents\\inlabru_workshop2025\\inlabru_workshop\\docs\\practicals\\GLM_ex.qmd"----

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


