mov.avg <- function(x,width=1000)
{
   ma <- numeric(length(x))
   for (i in 1:width)
      ma[i] <- mean(x[1:i])
   for (i in (width+1):length(x)) 
      ma[i] <- mean(x[(i-width):i])
   
   return(ma)
}

mov.stats <- function(x,width=1000)
{
   mw.stats <- matrix(NA,nrow=length(x),ncol=2)
   for (i in 1:width){
      mw.stats[i,1] <- mean(x[1:i])
      mw.stats[i,2] <- sd(x[1:i])
   }
   for (i in (width+1):length(x)){ 
      mw.stats[i,1] <- mean(x[(i-width):i])
      mw.stats[i,2] <- sd(x[(i-width):i])
   }
   
   return(mw.stats)
}

delta <- function(x)
{
   d <- numeric(length(x)-1)
   for (i in 2:length(x))
      d[i-1] <- x[i] - x[i-1]
   return(d)
}


div1 <- read.table("Diversity21.csv",header=T,sep=";")

nGen <- div1$Step/40000
nspec1000 <- mov.stats(div1$NSpec)
shannon1000 <- mov.stats(div1$Shannon)
pie1000 <- mov.stats(div1$NSpec)

par(mfrow=c(2,2))
plot(NSpec~nGen,data=div1,type="l",log="")
plot(Shannon~nGen,data=div1,type="l",log="")
lines(nGen,shannon1000[,1],col=2)
plot(PIE~nGen,data=div1,type="l",log="")

cor(div1[29000:30000,c(6,7)])

plot(nGen,shannon1000[,2],col=2,type="l")


# nspec100 <- mov.avg(div1$NSpec,100)
# nspec1000 <- mov.avg(div1$NSpec,1000)
# lines(nGen,nspec100,col=2)
# lines(nGen,nspec1000,col=3)


plot(Shannon~nGen,data=div1,type="l")
shannon100 <- mov.avg(div1$Shannon,100)
shannon1000 <- mov.avg(div1$Shannon,1000)

lines(nGen,shannon100,col=2)
lines(nGen,shannon1000,col=3)

ma10 <- mov.avg(div1$NSpec,width=10)
lines(nGen,ma10,col=2)
d10 <- delta(ma10)
abline(v=min(which(d10>0))*10)

acf(div1$NSpec[20:101],lag.max=10)



