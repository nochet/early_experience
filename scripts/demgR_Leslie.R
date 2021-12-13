###################################################
### chunk number 1: 
###################################################
library("demogR")
options(digits = 3, scipen = 1)  # biases toward fixed notation
data("goodman")
names(goodman)
ven <- with(goodman, life.table(x=age, nKx=ven.nKx, nDx=ven.nDx))
mad <- with(goodman, life.table(x=age, nKx=mad.nKx, nDx=mad.nDx))
usa <- with(goodman, life.table(x=age, nKx=usa.nKx, nDx=usa.nDx))
ven


###################################################
### chunk number 2: 
###################################################
data(thar)
thar.lt <- with(thar, life.table(x=age, nDx=deaths, nKx=count, type="cohort", iwidth=1, width12=c(1,1)))
thar.lt


###################################################
### chunk number 3: 
###################################################
ven.mx <- with(goodman, ven.bx/ven.nKx)
mad.mx <- with(goodman, mad.bx/mad.nKx)
A <- leslie.matrix(lx=ven$nLx,mx=ven.mx)
B <- leslie.matrix(lx=mad$nLx, mx=mad.mx)
A


###################################################
### chunk number 4: 
###################################################
Px <- c(0.77, 0.95, 0.98, 0.97)
Fx <- c(0,0,1,1.2,1)
L <- odiag(Px,-1)
L[1,] <- Fx
L


###################################################
### chunk number 5: 
###################################################
usa <- with(goodman, life.table(x=age, nKx=usa.nKx, nDx=usa.nDx))
usa.mx <- goodman$usa.bx/goodman$usa.nKx
C <- leslie.matrix(lx=usa$nLx,mx=usa.mx)
C


###################################################
### chunk number 6: 
###################################################
no <- goodman$usa.nKx[3:11]
no <- c(sum(goodman$usa.nKx[1:2]),no)
tmax <- 20
N <-  project.leslie(A=C,no=no, tmax=tmax)


###################################################
### chunk number 7: 
###################################################
cols <- rgb(0,(10:1)/10,(1:10)/10)
plot(5*(0:20),N[1,]/100000,
     type="l",
     xlab="Years", 
     ylab="Population Size (x100,000)",
     ylim=c(16,175), col=cols[1])
for(i in 2:10) lines(5*(0:20),N[i,]/100000, col=cols[i])


###################################################
### chunk number 8: 
###################################################
N.tot <- project.leslie(A=C,no=no, tmax=tmax, pop.sum=TRUE)
N.tot


###################################################
### chunk number 9: 
###################################################
Aea <- eigen.analysis(A)
Bea <- eigen.analysis(B)
Cea <- eigen.analysis(C)
Aea


###################################################
### chunk number 10: 
###################################################
calc.ro(A)
gen.time(A)


###################################################
### chunk number 11: 
###################################################
par(mfrow=c(2,3))
plot(Aea$sensitivities)
title("(a)")
plot(Bea$sensitivities)
title("(b)")
plot(Cea$sensitivities)
title("(c)")
plot(Aea$elasticities)
title("(d)")
plot(Bea$elasticities)
title("(e)")
plot(Cea$elasticities)
title("(f)")


###################################################
### chunk number 12: 
###################################################
secder(A,k=2,l=1)


###################################################
### chunk number 13: 
###################################################
elassens(A,2,1)


###################################################
### chunk number 14: 
###################################################
age <- seq(0,45,by=5)
SD <- fullsecder(B)
SD
cols <- rgb(0,(1:10)/10,(10:1)/10)
plot(age,SD[11,11:20],type="l", col=cols[1], xlab="Age", ylab="secder")
for(i in 12:20) lines(age,SD[i,11:20], col=cols[i-11])


###################################################
### chunk number 15: 
###################################################
ven.le <- loop.elas(A, draw.plot=FALSE)
mad.le <- loop.elas(B, draw.plot=FALSE)
usa.le <- loop.elas(C, draw.plot=FALSE)


###################################################
### chunk number 16: 
###################################################
par(mfrow=c(2,3))
barplot(ven.le, names.arg=5*1:11, horiz=TRUE, beside=TRUE, xlim=c(0,0.3), xlab="Loop Elasticity", ylab="Age")
barplot(mad.le, names.arg=5*1:11, horiz=TRUE, beside=TRUE, xlim=c(0,0.3), xlab="Loop Elasticity", ylab="Age")
barplot(usa.le, names.arg=5*1:10, horiz=TRUE, beside=TRUE, xlim=c(0,0.3), xlab="Loop Elasticity", ylab="Age")


###################################################
### chunk number 17: 
###################################################
cdw <- cdmltw() # defaults to female
names(cdw)
plot(cdw$age,cdw$lx[5,], type="l", lty=3, col="red",xlab="Age", 
ylab="Survivorship")
for(i in 6:14) lines(cdw$age,cdw$lx[i,], lty=3, col="red")
ache.nKx <- c(292, 1045, 827, 966, 817, 664, 502, 376, 327, 277, 199,
               180, 125, 83, 43,13, 2)
ache.nDx <- c(34, 46, 17, 10, 6, 8, 1, 6, 1, 4, 5, 3, 3, 3, 2, 1, 1)
ache.age <- c(0,1,seq(5,75,by=5))
ache.lt <- life.table(x=ache.age,nKx=ache.nKx,nDx=ache.nDx, type="cd")
kung.lx <- c(1.0000, 0.7900, 0.6715, 0.6045, 0.5925, 0.5925, 0.5510,
0.5400, 0.5080, 0.4475, 0.4045, 0.3765, 0.3350, 0.3250, 0.2275,
0.1595, 0.0990, 0.0390, 0.0200)
kung.age <- c(0,1,seq(5,85,by=5))
lines(ache.lt$x, ache.lt$lx, lwd=2)
lines(kung.age,kung.lx,lwd=2,col="grey")
legend(72,1,c("Ache","!Kung"), lwd=2, lty=1, col=c("black","grey"))


