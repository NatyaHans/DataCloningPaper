
library(diversitree)
library(ape)
set.seed(123)

# Read tree
phy<-read.tree("/blue/burleigh/nhans/GitHubNH/DataCloning/BiSSE/Data/whales_ultra.tre")
data<-read.csv("/blue/burleigh/nhans/GitHubNH/DataCloning/BiSSE/Data/bodysize_cetaceans_binarydata.csv",header=TRUE)

n.row<-1 # MLE from different models
mle.mat<-matrix(NA, nrow=n.row, ncol=6, dimnames=list(as.character(c("Bisse")),c("lambda0","lambda1","mu0","mu1","q01","q10" )))
n.step<-15000 # MCMC simulation steps


#checks
if(is.null(phy)){
    print("Null tree")
    } else{ print("Loaded the tree")}

if(!is.ultrametric(phy)){
    print("Non-ultrametric tree")
    } else{print("Ultrametric tree")}

if(!is.binary(phy)){
    print("Non binary tree")
    } else{print("Binary tree")}

#...........................

n.taxa<-length(phy$tip.label)


if(length(data[,1])!= n.taxa){
    print("Species in table do not match in number with species in the phylo")
    } else{("Yay we passed the checks")}

### Assiging tip states to phylogeny
phy$tip.state<-rep(NA,n.taxa)#empty tip state vector
counter<-0
for(i in 1:n.taxa){
    for (j in 1:n.taxa){
        if(phy$tip.label[i]==as.character(data[,1][j])){
            phy$tip.state[i]=data[,2][j]
            names(phy$tip.state)[i]=phy$tip.label[i]
            }   
        else{
           counter<-counter+1
           # print("Tip and species mismatch")  
            }
        }
}
## Checks
counter
save<-(n.taxa*n.taxa)-n.taxa ## number of times the above loop stays in the else side
if( counter!=save  || length(phy$tip.state)!= n.taxa){
    print("there are probably missing tip states")
    }else{
        print("everything checks out")}



# Fitting Bisse ML model
print("Running Bisse")
mytimeBisse <- format(Sys.time(), "%b_%d_%y")
state<-phy$tip.state
#myfile <- file.path(setwd("/scratch/lfs/nhans/00-DataCloning/03-BiasCorrection"), paste0("what",n.tree, "_", k , ".csv"))

lik.bisse<-make.bisse(phy,state)
p<-starting.point.bisse(phy)
fit.phy<-find.mle(lik.bisse, p)
mle.mat[1,]<-fit.phy$par

### Now Bayesian MCMC model fitting       
#k.vec<-c(2,4,8,16,32,64)
#k.vec<-c(2,4,8,16,32,64,128,256,512)

k.vec<-2
## Setting Priors
prior<-make.prior.exponential(1/ (2* (p[1] -p[3])))
tmp<- mcmc(lik.bisse, fit.phy$par, nsteps=n.step, prior=prior,lower=0, w=rep(1, 6), print.every=0)
w.phy<-diff(sapply(tmp[2:7], range)) #this is used in posterior info

# Data Cloning
i.vec<-vector(mode="numeric", length=n.row) 
# Note: This is an empty vector for MLE simulation in the other codes
# But we are using real phylogeny here so this is 1,since we are estimating 
# the ML on the real phylogeny only once

#k.vec<-c(2,4,8,16,32,64,128,256,512,1024,2048,4096) # Clone vector
tot<-length(i.vec)*length(k.vec)

#results.array <- array(0,c(n.step,6,length(i.vec))) # Empty array for Posteriors of each MLE
### This is only required if we have more than one MLE simulation
# SO used a matrix instead
#results.mcmc.mat<-matrix(NA,nrow=n.step,ncol=6,dimnames=list(as.character(1:n.step),c("lambda0","lambda1","mu0","mu1","q01","q10" )))
#results.clone<-matrix(NA,nrow=n.step,ncol=6,dimnames=list(as.character(1:n.step),c("lambda0","lambda1","mu0","mu1","q01","q10" )))

results.clone<- array(0,c(n.step,6,tot)) # Empty array for Posteriors of each MLE
## This empty array is similar to the one above
name.vec<-vector(mode="character", length=tot)  # Empty vector for MCMC simulation that work



#print("Starting MCMC")
# Posterior simulations for cetaceans
#samples<-mcmc(lik.bisse, fit.phy$par, nsteps=n.step, prior=prior,lower=0, w=w.phy, print.every=0)
#message.NH5<- paste("Saving results:", "whatever" , sep= "")
#print(message.NH5)
#results.mcmc.mat[,1:6]<- as.matrix(samples[,2:7])

### Data Cloning here
source("/blue/burleigh/nhans/GitHubNH/DataCloning/BiSSE/Scripts/MCMC_clone.R")
counter2<-0
for (k in k.vec){

    counter2<-counter2+1
    message.NH6<- paste("saving clones", k , sep= "")
    print(message.NH6)
    samples.clone<-mcmc(lik.bisse, fit.phy$par, nsteps=n.step, prior=prior,lower=0, w=w.phy, print.every=0, k.clone=k)
    results.clone[,1:6,counter2]<- as.matrix(samples.clone)[,2:7]
    name.vec[counter2]<-(paste("clone",k, sep=""))
}


#Naming the cloning data
dimnames(results.clone)[[3]]<-name.vec
dimnames(results.clone)[[2]]<-c("lambda0","lambda1","mu0","mu1","q01","q10" )
dimnames(results.clone)[[1]]<-c(as.character(1:n.step))
#write.table(results.clone,file = paste("cetaceans_clone_table_", mytimeBisse ,".txt",sep = ""))
write.csv(results.clone,file = paste("cetaceans_clone_",k.vec,"_", mytimeBisse ,".csv",sep = ""))

# Added for scaling eigen ( to get max eigen from without clones original posteriors)
#dimnames(results.mcmc.mat)[[3]]<-name.vec
#dimnames(results.mcmc.mat)[[2]]<-c("lambda0","lambda1","mu0","mu1","q01","q10" )
#dimnames(results.mcmc.mat)[[1]]<-c(as.character(1:n.step))
#write.csv(results.mcmc.mat,file = paste("cetaceans_posterior_table_", mytimeBisse ,".csv",sep = ""))




