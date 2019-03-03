############FUNCTIONS##########

source("ee_PreProcess_lifespan_functions.R")


#check for errors

library(tidyverse)

eelife <- read.table('../data/processed/ee_lifespan_only.txt',
                     sep="\t",header=TRUE,stringsAsFactors=FALSE) %>%
  select(-age)%>%
  rename(ID = id)


#calculate age (day)
eelife$setDate <- as.Date(eelife$setDate , "%m/%d/%y")
eelife$flipDate <- as.Date(eelife$flipDate , "%m/%d/%y")
eelife_age <- (eelife$flipDate - eelife$setDate)
eelife <- cbind(eelife,eelife_age)
eelife$age<-as.numeric(eelife$eelife_age)


# Create larval_adult diet column & compute age
eelife <- eelife  %>% 
  unite(larv_adult, larvalTreat, adultTreat, sep = "_", remove = FALSE)

eelife$age <- as.numeric(eelife$age)

eelife[320,]
hist(eelife$age)


#Find dupicated samples
duplnames <- duplicate.ages(unique(eelife$ID), eelife[,c('ID','age')])
duplnames

#eelife <- eelife[complete.cases(eelife[ , "age"]),]

#Age Check- skipped data entry (rows)

#Look for missing rows - Age gap >3 days
check.age <- age.check(unique(eelife$ID), eelife[,c('ID','age')])
dd<-eelife[,c('ID','age')]
check.age
write.table(check.age,file='../data/processed/omitedages.txt', 
            sep='\t',row.names = FALSE)


#change NewAge
eelife$NewAge<-eelife$age+2

#duplicate NewAges for a replicate

#get unique replicate ids
Uids<-unique(eelife$ID)

dups.all<-eelife[0,]
#loop through, find if duplicated, print ID and NewAge
for(ii in Uids)
{
  ID.dat<-subset(eelife, ID ==ii)
  if(anyDuplicated(ID.dat$NewAge)!=0)
  {
    dups.all<-rbind(dups.all,ID.dat[c(which(duplicated(ID.dat$NewAge)),which(duplicated(ID.dat$NewAge,fromLast=TRUE))),])
  }
  
}

dups.all<-dups.all[order(dups.all$ID,dups.all$NewAge),]
#write.table(dups.all, "duplicatedNewAges.txt",sep="\t",row.names=FALSE)



# construct NstatF and NstartM from the data
# For each ID, sum(deadF), sum(deadM), sum(cens)
StartNum <- eelife %>% 
  group_by(ID) %>% 
  summarise (sum_deadF = sum(deadF, n=n()),
             sum_deadM = sum(deadM, n=n()),
             sum_cens = sum(cens), n=n())

Start <- StartNum %>% 
  group_by(ID) %>% 
  summarise(Nstart = sum_deadF+sum_deadM+sum_cens)%>%
  as.data.frame()

# Fill-in NstarF with observed number dead
eelife$Nstart <- Start[match(eelife$ID, Start$ID),2] #2 selects the second column of F

eelife <- eelife %>%
  select(-eelife_age) %>%
  data.frame()

N.start<-616


#eelife must have columns ID, NewAge, Dead, Censored, Carried
proc.data<-Manip.Survival(eelife,N.start)
newdata<-proc.data$dat
nevents_by_ID<-proc.data$nevents
subset(nevents_by_ID, nevents>616)

newdata[1,]


save(newdata, file="newdata.rda") #processed original data
load("newdata.rda")

