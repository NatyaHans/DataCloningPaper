# This is the correct version of fitting Dipsidae
# Also made into an R script to run on cluster


library(diversitree)
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

if(!is.binary.tree(phy)){
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


## Setting Priors
prior<-make.prior.exponential(1/ (2* (p[1] -p[3])))
tmp<- mcmc(lik.bisse, fit.phy$par, nsteps=n.step, prior=prior,lower=0, w=rep(1, 6), print.every=0)
w.phy<-diff(sapply(tmp[2:7], range)) #this is used in posterior info

i.vec<-vector(mode="numeric", length=n.row) 
# Note: This is an empty vector for MLE simulation in the other codes
# But we are using real phylogeny here so this is 1,since we are estimating 
# the ML on the real phylogeny only once


#results.array <- array(0,c(n.step,6,length(i.vec))) # Empty array for Posteriors of each MLE
### This is only required if we have more than one MLE simulation
# SO used a matrix instead
results.mcmc.mat<-matrix(NA,nrow=n.step,ncol=6,dimnames=list(as.character(1:n.step),c("lambda0","lambda1","mu0","mu1","q01","q10" )))
#results.clone<-matrix(NA,nrow=n.step,ncol=6,dimnames=list(as.character(1:n.step),c("lambda0","lambda1","mu0","mu1","q01","q10" )))

#results.clone<- array(0,c(n.step,6,tot)) # Empty array for Posteriors of each MLE
## This empty array is similar to the one above
#name.vec<-vector(mode="character", length=tot)  # Empty vector for MCMC simulation that work



print("Starting MCMC")
# Posterior simulations for cetaceans
samples<-mcmc(lik.bisse, fit.phy$par, nsteps=n.step, prior=prior,lower=0, w=w.phy, print.every=0)
message.NH5<- paste("Saving results:", "No clone" , sep= "")
print(message.NH5)
results.mcmc.mat[,1:6]<- as.matrix(samples[,2:7])


# Added for scaling eigen ( to get max eigen from without clones original posteriors)
colnames(results.mcmc.mat)<-c("lambda0","lambda1","mu0","mu1","q01","q10" )
#dimnames(results.mcmc.mat)[[1]]<-c(as.character(1:n.step))
write.csv(results.mcmc.mat,file = paste("cetaceans_posterior_table_", mytimeBisse ,".csv",sep = ""))




