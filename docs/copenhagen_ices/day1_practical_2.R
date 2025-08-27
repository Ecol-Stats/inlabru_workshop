## ----child="C:/Users/admin/Documents/inlabru_workshop/inlabru_workshop/docs/practicals\\LM_ex_priors.qmd"------------------------

## --------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false

library(webexercises)



## --------------------------------------------------------------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     


## --------------------------------------------------------------------------------------------------------------------------------
#| code-fold: show
beta = c(2,0.5)
sd_error = 0.1

n = 100
x = rnorm(n)
y = beta[1] + beta[2] * x + rnorm(n, sd = sd_error)

df = data.frame(y = y, x = x)  



## --------------------------------------------------------------------------------------------------------------------------------
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


## --------------------------------------------------------------------------------------------------------------------------------
#| eval: true
#| purl: true 

inla.priors.used(fit.lm)






## --------------------------------------------------------------------------------------------------------------------------------

# First we define the logGamma (0.01,0.01) prior 

prec.tau <- list(prec = list(prior = "loggamma",   # prior name
                             param = c(0.01, 0.01))) # prior values

lik2 =  bru_obs(formula = y ~.,
                family = "gaussian",
                data = df,
                control.family = list(hyper = prec.tau))

fit.lm2 = bru(cmp2, lik2) 



## --------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 4
#| fig-align: center
#| fig-height: 4
#| code-fold: show

plot(fit.lm, "beta_0")




## --------------------------------------------------------------------------------------------------------------------------------
res_samples <- predict(
  fit.lm,         # the fitted model
  df,             # the original data set
  ~ data.frame(   
    res = y-(beta_0 + beta_1)  # compute the residuals
  ),
  n.samples = 1000   # draw 1000 samples
)



## --------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 6
#| fig-height: 4
#| fig-align: center
#| fig-cap: "Bayesian residual plots: the left panel is the residual index plot; the right panel is the plot of the residual versus the covariate x"
#| code-fold: true
#| code-summary: "Residuals checks for Linear Model"

ggplot(res_samples,aes(y=mean,x=1:100))+geom_point() +
ggplot(res_samples,aes(y=mean,x=x))+geom_point()



## --------------------------------------------------------------------------------------------------------------------------------
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



## ----child="C:/Users/admin/Documents/inlabru_workshop/inlabru_workshop/docs/practicals\\LMM_fish_example.qmd"--------------------

## --------------------------------------------------------------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     



## --------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false
#| eval: false

# library(FSAdata)
# df = PygmyWFBC
# # subset for 2001 and remove nets with zero observations
# df = df %>% filter(year==2001 & !is.na(sex) & net_no %in% c(1:5,7:8)) %>% dplyr::select(net_no,wt,tl,sex)
# write.csv(df,file = "datasets/PygmyWFBC.csv",row.names = F)
# 


## --------------------------------------------------------------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
#| warning: false
#| message: false

PygmyWFBC <- read.csv("datasets/PygmyWFBC.csv")

ggplot(PygmyWFBC, aes(x = factor(net_no), y = wt,fill = sex)) + 
  geom_boxplot() + 
  labs(y="Weight (g)",x = "Net no.")



## --------------------------------------------------------------------------------------------------------------------------------
PygmyWFBC$sex_M <- ifelse(PygmyWFBC$sex=="F",0,1)


## --------------------------------------------------------------------------------------------------------------------------------
cmp =  ~ -1 + sex_M +  beta_0(1)  + beta_1(tl, model = "linear") +   net_eff(net_no, model = "iid")

lik =  bru_obs(formula = wt ~ .,
            family = "gaussian",
            data = PygmyWFBC)

fit = bru(cmp, lik)

summary(fit)




## --------------------------------------------------------------------------------------------------------------------------------
plot(fit,"net_eff")



## --------------------------------------------------------------------------------------------------------------------------------

sampvars <- 1/inla.hyperpar.sample(1000,fit,improve.marginals = T)

colnames(sampvars) <- c("Error variance","Between-net Variance")

apply(sampvars,2,
      function(x) c("mean"=mean(x),
                    "std.dev" = sd(x),
                    quantile(x,c(0.025,0.5,0.975))))



## ----child="C:/Users/admin/Documents/inlabru_workshop/inlabru_workshop/docs/practicals\\GLM_scores_ex.qmd"-----------------------

## --------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false

library(webexercises)



## --------------------------------------------------------------------------------------------------------------------------------
#| warning: false
#| message: false

library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     



## --------------------------------------------------------------------------------------------------------------------------------
crabs <- read.csv("datasets/crabs.csv")

# conditional means and variances
crabs %>%
  summarise( Mean = mean(satell ),
             Variance = var(satell),
                     .by = color)


## --------------------------------------------------------------------------------------------------------------------------------

crabs_df = model.matrix( ~  color , crabs) %>%
  as.data.frame() %>%
  select(-1) %>%        # drop intercept
  bind_cols(crabs) %>%  # append to original data
  select(-color)        # remove original color categorical variable



## --------------------------------------------------------------------------------------------------------------------------------

cmp =  ~ -1 + beta0(1) +  colordarker +
       colorlight + colormedium +
       w(weight, model = "linear")

lik =  bru_obs(formula = satell ~.,
            family = "poisson",
            data = crabs_df)

fit_pois = bru(cmp, lik)

summary(fit_pois)



## --------------------------------------------------------------------------------------------------------------------------------
bru_options_set(control.compute = list(cpo = TRUE))

fit_pois = bru(cmp, lik)




## --------------------------------------------------------------------------------------------------------------------------------
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



## ----child="C:/Users/admin/Documents/inlabru_workshop/inlabru_workshop/docs/practicals\\HGAMM_ex.qmd"----------------------------

## --------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false

library(webexercises)



## --------------------------------------------------------------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     



## --------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| fig-align: center
#| fig-width: 6
#| fig-height: 6
#| fig-cap: "Source–depth proﬁles per month. Each line represents a station."
icit <- read.csv("datasets/ISIT.csv")

icit$Month <- as.factor(icit$Month)
levels(icit$Month) <- month.abb[unique(icit$Month)]
icit$Month_id <- as.numeric(icit$Month)

ggplot(icit,aes(x=SampleDepth,y= Sources,group=as.factor(Station),colour=as.factor(Station)))+geom_line()+facet_wrap(~Month)+theme(legend.position = "none")


## --------------------------------------------------------------------------------------------------------------------------------
#| eval: false

# icit <- read.csv("datasets/ISIT.csv")
# 
# icit$Month <- as.factor(icit$Month)
# levels(icit$Month) <- month.abb[unique(icit$Month)]
# 
# ggplot(icit,aes(x=SampleDepth,y= Sources,
#                 group=as.factor(Station),
#                 colour=as.factor(Station)))+
#   geom_line()+
#   facet_wrap(~Month)+
#   theme(legend.position = "none")


## --------------------------------------------------------------------------------------------------------------------------------

icit$Month_id <- as.numeric(icit$Month) # numeric index for the i-th month

cmp_g =  ~ -1+ beta_0(1) + 
  smooth_g(SampleDepth, model = "rw1") + 
  month_reff(Month_id, model = "iid") 

lik =  bru_obs(formula = Sources ~.,
               family = "gaussian",
               data = icit)

fit_g = bru(cmp_g, lik)

summary(fit_g)



## --------------------------------------------------------------------------------------------------------------------------------
#| code-summary: "Global smoother marginal effect"
#| code-fold: true
#| fig-align: center
#| fig-width: 4
#| fig-height: 4

data.frame(fit_g$summary.random$smooth_g) %>% 
  ggplot() + 
  geom_ribbon(aes(ID,ymin = X0.025quant, ymax= X0.975quant), alpha = 0.5) + 
  geom_line(aes(ID,mean)) + 
  xlab("covariate") + ylab("smooth effect")


## --------------------------------------------------------------------------------------------------------------------------------
icit$depth_grouped <- inla.group(icit$SampleDepth,n=50)





## --------------------------------------------------------------------------------------------------------------------------------

cmp_gs =  ~ -1+ beta_0(1) +
  smooth_g(SampleDepth, model = "rw1") + 
  month_reff(Month_id, model = "iid")+
  smooth_loc(SampleDepth, model = "rw1", group = Month_id)



## --------------------------------------------------------------------------------------------------------------------------------

fit_gs = bru(cmp_gs, lik) 



## --------------------------------------------------------------------------------------------------------------------------------
pred_gs = predict(fit_gs, icit, ~ (beta_0 + smooth_g+month_reff+smooth_loc))



## --------------------------------------------------------------------------------------------------------------------------------
#| code-summary: "Global + group smoother predictions"
#| code-fold: true

ggplot(pred_gs,aes(y=mean,x=SampleDepth))+
  geom_ribbon(aes(SampleDepth,ymin = q0.025, ymax= q0.975), alpha = 0.5,fill="tomato") +
  geom_line()+
  geom_point(aes(x=SampleDepth,y=Sources ),alpha=0.25,col="grey40")+
  facet_wrap(~Month)




