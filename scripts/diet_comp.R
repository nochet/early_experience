
library(tidyverse)
library(data.table)

dcomp <- read.csv("/Users/ngomae/MyGithub/early_experience/processed/diet_composition.csv") 

dcomp <-as.data.frame(t(dcomp))
names(dcomp) <- dcomp[1,] 
dcomp <- dcomp[2:5,]


