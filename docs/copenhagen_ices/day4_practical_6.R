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


