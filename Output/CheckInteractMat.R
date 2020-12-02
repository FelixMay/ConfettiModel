mat <- as.matrix(read.table("InteractMat_sim1_rep1.csv",header=F,sep=";"))
hist1 <-hist(mat)

 hist(diag(mat))
mean(diag(mat))
sd(diag(mat))
dim(mat)

spec.dat <- read.table("SpeciesOut_sim1_rep1.txt",header=T,sep="\t",row.names=NULL)
dim(spec.dat)
head(spec.dat)

plot(spec.dat$trait, spec.dat$LocalAbund)

d1 <- as.matrix(dist(spec.dat$trait))
plot(d1, mat)

hist(spec.dat$trait)

plot(diag(mat),spec.dat$CNDD)
abline(0,1)

plot(Kcon50 ~ Kcon10, data=spec.dat,log="xy")
plot(Kcon50 ~ LocalAbund, data=spec.dat,log="xy")

plot(Kcon10 ~ CNDD, data=spec.dat,log="y")

plot(Kcon50 ~ muDisp,data=spec.dat,log="y")
cor.test(spec.dat$Kcon10[spec.dat$Kcon10>0],spec.dat$muDisp[spec.dat$Kcon10>0],method="spearman")

summary(spec.dat)
sum(spec.dat$LocalAbund)

abund <- read.table("Abund1.csv",header=F,sep=";")
head(abund)
rowSums(abund)
