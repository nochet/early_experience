
source("PreProcess_lifespan_functions.R")


library(tidyverse)

eelife <- read.table("../processed/ee_lifespan_only.txt",
                     sep="\t",header=TRUE,stringsAsFactors=FALSE) %>%
  select(-age)

#calculate age (day)
eelife$setDate <- as.Date(eelife$setDate , "%m/%d/%y")
eelife$flipDate <- as.Date(eelife$flipDate , "%m/%d/%y")
eelife_age <- (eelife$flipDate - eelife$setDate)
eelife <- cbind(eelife,eelife_age)
eelife$age<-as.numeric(eelife$eelife_age)
#all.equal(eelife$age, eelife$days)
#wh.diff<-eelife$age-eelife$days

# Create larval_adult diet column & compute age
eelife <- eelife  %>% 
  unite(larv_adult, larvalTreat, adultTreat, sep = "_", remove = FALSE)

eelife$age <- as.numeric(eelife$age)

eelife[320,]
hist(eelife$age)

#Make unique ids
#eelife$id <- paste(eelife$fID,'_',eelife$treat,sep='')

#Find dupicated samples
duplnames <- duplicate.ages(unique(eelife$id), eelife[,c('id','age')])
duplnames

#eelife <- eelife[complete.cases(eelife[ , "age"]),]

#Age Check- skipped data entry (rows)

#Look for missing rows - Age gap >3 days
check.age <- age.check(unique(eelife$id), eelife[,c('id','age')])
dd<-eelife[,c('id','age')]
check.age
write.table(check.age,file='../processed/omitedages.txt', 
            sep='\t',row.names = FALSE)

#change NewAge
eelife$NewAge<-eelife$age+2

#Check if some letters in fID are in lower case (i.e. expect all names in upper case)
eelife %in% letters

#write out cleaned data
write.table(eelife,file='../processed/eelife_correctedData.txt', 
            sep='\t',row.names = FALSE)

################################


eelife<-read.table("../processed/eelife_correctedData.txt",
                   sep="\t",stringsAsFactors=FALSE,header=TRUE)

# Create 'carriedF' and 'carriedM' for function to run, fill with 0s
eelife$carriedF <- rep(0,nrow(eelife)) 
eelife$carriedM <- rep(0,nrow(eelife))

#separate all events into rows
#females
Flife.dat<-eelife[,c('setDate','flipDate','age','NewAge','id','larv_adult',
                     'larvalTreat','adultTreat','deadF','carriedF')]
colnames(Flife.dat)[which(colnames(Flife.dat)=='deadF')]<-'Dead'
colnames(Flife.dat)[which(colnames(Flife.dat)=='carriedF')]<-'Carried'
#Flife.dat<-Flife.dat[-which(is.na(Flife.dat$Dead)),]
Fevent<-Manip.Survival(Flife.dat)

#males
Mlife.dat<-eelife[,c('setDate','flipDate','age','NewAge','id','larv_adult',
                     'larvalTreat','adultTreat','deadM','carriedM')]
colnames(Mlife.dat)[which(colnames(Mlife.dat)=='deadM')]<-'Dead'
colnames(Mlife.dat)[which(colnames(Mlife.dat)=='carriedM')]<-'Carried'
#Mlife.dat<-Mlife.dat[-which(is.na(Mlife.dat$Dead)),]
Mevent<-Manip.Survival(Mlife.dat)

#censored events
life.dat<-eelife[,c('setDate','flipDate','NewAge','id','cens')]
colnames(life.dat)[which(colnames(life.dat)=='cens')]<-'Censored'
#life.dat<-life.dat[-which(is.na(life.dat$Censored)),]
Cevent<-Cen.events(life.dat)

####

# construct NstatF and NstartM from the data
# For each id, sum(deadF), sum(deadM), sum(cens)
StartNum <- eelife %>% 
  group_by(id) %>% 
  summarise (sum_deadF = sum(deadF, n=n()),
             sum_deadM = sum(deadM, n=n()),
             sum_cens = sum(cens), n=n())

F <- StartNum %>%
  select(id, sum_deadF)%>%
  as.data.frame()

M <- StartNum %>%
  select(id, sum_deadM) %>%
  as.data.frame()

# Fill-in NstarF with observed number dead
eelife$NstartF <- F[match(eelife$id, F$id),2] #2 selects the second column of F
eelife$NstartM <- M[match(eelife$id, F$id),2]

eelife <- eelife %>%
  select(-eelife_age) %>%
  data.frame()
#colnames(eelife)[colnames(eelife)=="NewAge"] <- "age"

totals <- CountEvents(eelife)
max(totals$NCensor)
which.max(totals$NCensor)
totals[53, ]   # known case of escapes
subset(eelife, id==totals[which.max(totals$NCensor),'id'] & cens>0)

################################
#Account for censored events
totals$miss.f <- totals$NstartF-totals$NdeadF
totals$miss.m <- totals$NstartM-totals$NdeadM  

for(jj in 1:nrow(Cevent)) 
{
  miss.f<-totals[totals$id==Cevent$id[jj],'miss.f']
  miss.m<-totals[totals$id==Cevent$id[jj],'miss.m']
  
  
  if(miss.f>miss.m)
    {
      Fcol<-Fevent[Fevent$id==Cevent$id[jj],-which(colnames(Fevent) %in% colnames(Cevent))][1,]
      NewFev<-cbind(Fcol,Cevent[jj,])
      NewFev<-NewFev[,colnames(Fevent)]
      Fevent<-rbind(Fevent,NewFev)
      totals$miss.f[totals$id==Cevent$id[jj]]<-totals$miss.f[totals$id==Cevent$id[jj]]-1
    }else{
      if(miss.m>miss.f)
      {
        Mcol<-Mevent[Mevent$id==Cevent$id[jj],-which(colnames(Mevent) %in% colnames(Cevent))][1,]
        NewMev<-cbind(Mcol,Cevent[jj,])
        NewMev<-NewMev[,colnames(Mevent)]
        Mevent<-rbind(Mevent,NewMev)
        totals$miss.m[totals$id==Cevent$id[jj]]<-totals$miss.m[totals$id==Cevent$id[jj]]-1
        
      }else{
        #randomly choose
        if(sample(c(1,2),1)==1)
        {
          Fcol<-Fevent[Fevent$id==Cevent$id[jj],-which(colnames(Fevent) %in% colnames(Cevent))][1,]
          NewFev<-cbind(Fcol,Cevent[jj,])
          NewFev<-NewFev[,colnames(Fevent)]
          Fevent<-rbind(Fevent,NewFev)
          totals$miss.f[totals$id==Cevent$id[jj]]<-totals$miss.f[totals$id==Cevent$id[jj]]-1
          
        }else{
          Mcol<-Mevent[Mevent$id==Cevent$id[jj],-which(colnames(Mevent) %in% colnames(Cevent))][1,]
          NewMev<-cbind(Mcol,Cevent[jj,])
          NewMev<-NewMev[,colnames(Mevent)]
          Mevent<-rbind(Mevent,NewMev)
          totals$miss.m[totals$id==Cevent$id[jj]]<-totals$miss.m[totals$id==Cevent$id[jj]]-1
          
          
        }
      }
    }
      
   
  
}

length(Mevent$status[Mevent$status==3])+length(Fevent$status[Fevent$status==3])


totals$miss.m

Fevent[1,]

write.table(Fevent, "../processed/Female_events_eelife.txt", row.names=FALSE, sep="\t")
write.table(Mevent, "../processed/Male_events_eelife.txt", row.names=FALSE, sep="\t")

## End of script

