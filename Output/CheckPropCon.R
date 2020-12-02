fr1 <- read.table("PropCon1.csv",header=T,sep=";")
plot(1:50,fr1[1,],type="l",ylim=c(0,max(fr1)))
for (i in 1:6) lines(1:50,fr1[i,],type="l",col=i)
#for (i in 4:6) lines(1:50,fr1[i,],type="l",col=3)


gr <- read.table("PCF1.csv",header=T,sep=";")
plot(1:50,gr[1,],type="l",ylim=c(0,1.2))
for (i in 1:6) lines(1:50,gr[i,],type="l",col=i)
abline(v=5)
for (i in 4:6) lines(1:50,gr[i,],type="l",col=i)


trees <- read.table("Trees_sim1_rep1.csv",header=T,sep=";")
rb.col <- rainbow(length(table(trees$SpecID)))
plot(Y~X,col=rb.col[SpecID],data=subset(trees,X<100 & Y<100))


rdist = sqrt(pi/2.0) * abs(rnorm(100000,0,20))
hist(rdist)
summary(rdist)

dd <- read.table("Test.txt")
dd <- as.numeric(dd[[1]])
hist(dd,freq=F)
lines(density(rdist),col=2)
mean(dd)

library(spatstat)

pp1 <- ppp(x=trees$X,y=trees$Y,window=owin(c(0,500),c(0,500)))
pcf1 <- pcf(pp1)
plot(pcf1)
lines(seq(1,100,by=1),gr[1,],col=3)

meanDisp <- 30
sdDisp <- 30

sigmaDisp = sqrt(log(1.0 + (sdDisp^2)/(meanDisp^2)))
muDisp = log(meanDisp) - 0.5 * sigmaDisp^2

x <- seq(0,100,by=0.001)
y <- dlnorm(x,meanlog=muDisp,sdlog=sigmaDisp)
plot(x,y,type="l")
