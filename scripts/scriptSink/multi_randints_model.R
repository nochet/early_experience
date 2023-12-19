
library(dplyr)
library(rstanarm)

eelife <- read.csv("/home/engoma/MyGithub/early_experience/processed/eeAlldat.csv")
eelife$status <- ifelse(eelife[,"status"] == 2, 1, 0)
eelife <- select(eelife, -X)

start_time <- Sys.time()

seed <- set.seed(717871)
chains <- 6
cores  <- 6
#iter   <- 1000

mod4 <- stan_surv(formula = Surv(age, status) ~ 
    larvalTreat + adultTreat + sex  + 
    (larvalTreat | cageID) + 
    (adultTreat | cageID) + 
    (sex | cageID),
  data    = eelife,
  basehaz = "ms", basehaz_ops = list(df = 10),
  chains  = chains, 
  cores   = cores, 
  seed    = seed,
  warmup = 15000,
  iter    = 5e5)

end_time <- Sys.time()
diff.t <- end_time-start_time

save(mod4, file="../processed/cageID_multi_randints.Rda")


