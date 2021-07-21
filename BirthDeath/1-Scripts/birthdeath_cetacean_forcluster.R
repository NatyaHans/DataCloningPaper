library(diversitree)
library(ape)
library(phytools)
set.seed(123)
mytimebd <- format(Sys.time(), "%b_%d_%y")

phy.name<-read.tree("../Data/whales_ultra.tre")
lik.phy<-make.bd(phy.name)
parML<-fit.bd(phy.name)

# without specifying priors
n.step<-15000
k.vec<- 2^seq(1,9,by=1)
samples.noprior<-array(0,c(n.step,2,length(k.vec)))# Empty array for Posteriors of each MLE

#----- No cloning
library(ape)
print("No cloning")
p<-starting.point.bd(phy.name)
prior.phy<-make.prior.exponential(1/ (2* (p[1] -p[2])))
samples.noclone<-mcmc(lik.phy, c(parML$b,parML$d), nsteps=n.step, prior=prior.phy,lower=0, w = c(0.1, 0.1), print.every=1000)
max.eigen.bd<-max(eigen(var(samples.noclone[2:3]))$values)
# this max.eigen is to be used to scale eigen mat
#------

source("/blue/burleigh/nhans/GitHubNH/DataCloning/BiSSE/Scripts/MCMC_clone.R")
#source("MCMC_clone.R")
counter<-0
for(k in k.vec){
counter<-counter+1
print(k)
assign(paste("samplesnoprior", k ,sep=""),mcmc(lik.phy, c(parML$b,parML$d), nsteps = n.step, lower = c(-Inf, -Inf), upper = c(Inf, Inf), w = c(0.1, 0.1), fail.value = -Inf, print.every = 1000,k.clone=k))
samples.noprior[,,counter]<- as.matrix(eval(parse(text=paste("samplesnoprior",  k, sep="")))[2:3])
}

# specifying exponential priors
#samples.withprior<-matrix(0,nrow=n.step, ncol=2)# Empty array for Posteriors of each MLE
samples.withprior<-array(0,c(n.step,2,length(k.vec)),dimnames=list(as.character(1:n.step),c("lambda","mu"),as.character(k.vec)))
p<-starting.point.bd(phy.name)
prior.phy<-make.prior.exponential(1/ (2* (p[1] -p[2])))
counter<-0
for(k in k.vec){
counter<-counter+1
print(k)
assign(paste("samples", k, sep=""),mcmc(lik.phy, c(parML$b,parML$d), nsteps=n.step, prior=prior.phy,lower=0, w = c(0.1, 0.1), print.every=1000, k.clone=k))
samples.withprior[,,counter]<- as.matrix(eval(parse(text=paste("samples", k,sep="")))[2:3])
}
write.csv(samples.noclone,file="posteriorbd_15000_Jul21.csv")
write.csv(samples.withprior,file="clonebd_15000_Jul21.csv")
