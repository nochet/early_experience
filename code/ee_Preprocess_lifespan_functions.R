

Manip.Survival<-function(lifedat,N.start)
#lifedat must have columns ID, NewAge, Dead, Censored, Carried
  {

  #get event indicies
D.events<-which(lifedat$Dead!=0)
Cen.events<-which(lifedat$Censored!=0)

#set up output
lifeInd<-lifedat[1,-which(colnames(lifedat) %in% c('Dead','Censored','Carried'))]
lifeInd$status<-0
lifeInd<-lifeInd[0,]

for (i in 1:length(D.events))
{
  ss<-lifedat[D.events[i],]
  NewAges<-sort(unique(lifedat[lifedat$ID==ss$ID,'NewAge']))
  if(ss$NewAge>min(NewAges))
  {
    #get data from previous NewAge
    #add NA if
    ss.prev<-subset(lifedat, ID==ss$ID & NewAge==NewAges[which(NewAges==ss$NewAge)-1])
    #print error if duplicated
    if(nrow(ss.prev)>1){stop("Duplicated NewAges-Run Data Check")}
    
  nDead<-ss$Dead-ss.prev$Carried
  
  }else{
    nDead<-ss$Dead
  }
  if(nDead>0)
  {
    dd<-lifedat[D.events[i],-which(colnames(lifedat) %in% c('Dead','Censored','Carried'))]
    #this is a way to replicate rows of a data frame
    #it is kind of like doing this dd<-dd[c(7,7,7,7,7),]
    dd<-dd[rep(row.names(dd),nDead),]
    
    #2 = dead
    dd$status<-rep(2,nDead)
    lifeInd<-rbind(lifeInd,dd)
  }

}


d.cen<-lifedat[Cen.events,]
d.cen<-d.cen[rep(row.names(d.cen),d.cen[,'Censored']),-which(colnames(d.cen) %in% c('Dead','Censored','Carried'))]
d.cen$status<-rep(3,nrow(d.cen))
  
lifeInd<-rbind(lifeInd,d.cen)

#still alive

ids<-unique(lifedat$ID)

err.ids<-data.frame("id"=character(length(ids)), "nevents"=numeric(length(ids)),stringsAsFactors=FALSE)

for(j in 1:length(ids))
{
  tt<-subset(lifedat, ID==ids[j])
  n.event<-length(lifeInd[lifeInd$ID==ids[j],'status'])
  
  err.ids[j,]<-c(ids[j],n.event)
  
  if(N.start<n.event)
  {
    N.total<-n.event
  }else{
    N.total<-N.start
  }
  if(n.event<N.start)
  {
    d.rest<-tt[,-which(colnames(tt) %in% c('Dead','Censored','Carried'))]
    d.rest<-d.rest[nrow(d.rest),]
    d.rest<-d.rest[rep(row.names(d.rest),N.total-n.event),]
    d.rest$status<-1
    lifeInd<-rbind(lifeInd,d.rest)

  }
}
return(list('dat'=lifeInd,'nevents'=err.ids))

}

#function to look for duplicated ids

duplicate.ages<-function(Uids,dd)
  #Uids is unique vial id
  #dd is 2 columns, id is in column 1, age is in column 2
{
  dups.all<-dd[0,]
  #loop through, find if duplicated, print ID and age
  for(ii in Uids)
  {
    id.dat<-dd[dd[,1]==ii,]
    if(anyDuplicated(id.dat[,2])!=0)
    {
      dups.all<-rbind(dups.all,id.dat[c(which(duplicated(id.dat[,2])),which(duplicated(id.dat[,2],fromLast=TRUE))),])
      
    }
  }
  return(dups.all)
}

#check missing rows - identify these by age gap > 3 days

age.check <-function(Uids,dd)
  #Uids is unique vial id
  #dd is 2 columns, id is in column 1, age is in column 2
{
  dups.all<-dd[0,]
  #loop through, find if duplicated, print ID and age
  for(ii in Uids)
  {
    id.dat<-dd[dd[,1]==ii,]
    id.dat <- id.dat[order(id.dat[,2]),]
    max.diff <- max(diff(id.dat[,2]))
    if(max.diff > 3)
    {
      ww <- which(diff(id.dat[,2])>3)
      ww <- sort(c(ww,ww+1))
      dups.all<-rbind(dups.all,id.dat[ww,])
      
    }
    cat(ii,"\n")
  }
  return(dups.all)
}


Cen.events<-function(lifedat)
  #lifedat must have columns id, NewAge, Censored
{
  
  #get event indicies
  Cen.events<-which(lifedat$Censored!=0)
  
  #set up output
  lifeInd<-lifedat[1,-which(colnames(lifedat) %in% c('Censored'))]
  lifeInd$status<-0
  lifeInd<-lifeInd[0,]
  
  d.cen<-lifedat[Cen.events,]
  d.cen<-d.cen[rep(row.names(d.cen),d.cen[,'Censored']),-which(colnames(d.cen) %in% c('Censored'))]
  d.cen$status<-rep(3,nrow(d.cen))
  
  lifeInd<-rbind(lifeInd,d.cen)
  return(lifeInd)
  
}
