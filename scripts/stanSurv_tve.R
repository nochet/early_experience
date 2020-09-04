# Load packages

library(survival)
library(rstanarm)
library(tidyverse)

eelife <- read.csv("../processed/eeAlldat.csv") %>%
  select(-X)


# Bayesian survival models

start_time <- Sys.time()

### M-splines is the default in stan_surv() 

# Baseline hazard function - assuming covariate effect of zero
mod1 <- stan_surv(formula = Surv(age, status) ~ tve(larvalTreat) + 
                    tve(adultTreat) + tve(sex), 
                  data    = eelife,
                  chains  = 8, 
                  cores   = 4, 
                  seed    = 414867,
                  iter    = 1e5)
print(mod1, digits = 3)
# 1000 iter for 30 min

# Test different parametric hazards
# i.e. fit several models, each with a different baseline hazard

mod1_exp      <- update(mod0, basehaz = "exp") 
mod1_weibull  <- update(mod0, basehaz = "weibull") 
mod1_gompertz <- update(mod0, basehaz = "gompertz") 
mod1_bspline  <- update(mod0, basehaz = "bs") 
mod1_mspline1 <- update(mod0, basehaz = "ms") 
mod1_mspline2 <- update(mod0, basehaz = "ms", 
                        basehaz_ops = list(df = 10))

save(mod0,mod1_exp,mod1_weibull,mod1_gompertz,mod1_bspline,
        mod1_mspline1,mod1_mspline2, file="../processed/BaseHaz_new.Rda")

# Compare fits using `loo` for `stansurv` objects

lcomp <- loo_compare(loo(mod1_exp),
                     loo(mod1_weibull),
                     loo(mod1_gompertz),
                     loo(mod1_bspline),
                     loo(mod1_mspline1),
                     loo(mod1_mspline2))

# Plot the baseline hazards with 95% posterior uncertainty limits
load("../processed/BaseHaz.Rda")

plotfun <- function(model, title) {
  plot(model, plotfun = "basehaz") +              
    coord_cartesian(ylim = c(0,0.4)) +            
    labs(title = title) +                        
    theme(plot.title = element_text(hjust = 0.5,size = 11),
          axis.title = element_text(size = 10)) + 
    theme_set(theme_cowplot()) +
    theme_half_open() 
}
p_exp      <- plotfun(mod1_exp,      title = "Exponential")
p_weibull  <- plotfun(mod1_weibull,  title = "Weibull")
p_gompertz <- plotfun(mod1_gompertz, title = "Gompertz")
p_bspline  <- plotfun(mod1_bspline,  title = "B-splines with df = 5")
p_mspline1 <- plotfun(mod1_mspline1, title = "M-splines with df = 5")
p_mspline2 <- plotfun(mod1_mspline2, title = "M-splines with df = 10")
cs <-bayesplot::bayesplot_grid(p_exp,
                               p_weibull,
                               p_gompertz,
                               p_bspline,
                               p_mspline1,
                               p_mspline2,
                               grid_args = list(ncol = 3))

ggplot2::ggsave(cs,filename="../plots/baseHaz_tve_fixedEffs_new.pdf",height = 6, width=10)


## Record execution time and multicore use
end_time <- Sys.time()
diff_time <- end_time - start_time

