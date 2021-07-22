#-----------------------------------------------------------#
#Script for simulating trees hidden state using diversitree #
#------------------------------------------------------------#
# loading packages
set.seed(123)
library(diversitree)
library(numDeriv)
suppressWarnings(library(hisse)) # warnings added
mytimeHisse <- format(Sys.time(), "%b_%d_%y")


# Read the data
tree<-read.tree("../Data/whales_ultra.tre")
data<-read.csv("../Data/bodysize_cetaceans_binarydata_forHisse.csv") # for some reason this works. whatever
n.tree<-"whales"
#Hidden state associated with both of the observed state?
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)
# transition rate matrix
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
trans.rates.nodual.allequal = trans.rates.nodual
trans.rates.nodual.allequal[!is.na(trans.rates.nodual.allequal) &
                              !trans.rates.nodual.allequal == 0] = 1
trans.rates.nodual.allequal
#---------------------------------------------------------------------#

# MCMC 
n.step<-15000 #MCMC simulation
n.save<-n.step/100
k.vec<-2 # Clone vector

#pp = hisse(tree, data, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

# MCMC no clone
lik.hisse <- makeHiSSELikelihood(phy = tree, data = data, hidden.states = TRUE)
lik.func <- lik.hisse$log.lik
temp.p<-starting.point.musse(tree,4)

# note that the above function spits out the starting point parameters in lambda1.......mu4....q12.....format
temp.prior <-make.prior.exponential(1/(temp.p[1]+temp.p[3]))
#samples.noclone<-mcmc(lik.func, temp.p, nsteps=n.step, prior=temp.prior,lower=0, w=rep(1,20), print.every=0, save.every=100,save.file=paste("posteriorsMCMC", n.tree ,"_",mytimeHisse,"_",n.step,".csv",sep = ""))

#count<-0
#Simulations for Posteriors-MCMC
#message.NH5<- paste("Saving results:", count , sep= "")
#print(message.NH5)
#write.csv(samples.noclone,file=paste("posteriors", n.tree ,"_",mytimeHisse,"_",n.step,".csv",sep = "")) 

# Saving clone results

results.clone<- array(0,c(n.save,length(temp.p),length(k.vec))) # Empty array for Posteriors of each MLE
name.vec<-vector(mode="character", length=length(k.vec))  # Empty vector for MCMC simulation that work


#2) Doing cloning  here 
counter<-0
source("/blue/burleigh/nhans/GitHubNH/DataCloning/BiSSE/Scripts/MCMC_clone.R")
for (k in k.vec){
   counter<-counter+1
   message.NH6<- paste("saving clones", k , sep= "")
   print(message.NH6)
   assign(paste("samples_clone_", k, sep=""), mcmc(lik.func, temp.p, nsteps=n.step, prior=temp.prior, lower=0, w=rep(1:20), print.every=0, k.clone=k,save.every=100,save.file=paste("clone_taxa_MCMC", n.tree ,"_",k,"_",mytimeHisse,"_",n.step,".csv",sep = "")))
   results.clone[,,counter]<- as.matrix(eval(parse(text=paste("samples_clone_", k, sep="")))[,2:21])
   name.vec[counter]<-(paste("samples_clone_", k, sep=""))
   }
print("writing results")
#parnames<-c("lambda1","lambda2", "lambda3", "lambda4", "mu1","mu2", "mu3", "mu4", "q12","q13","q14","q21","q23","q24","q31","q32","q34","q41","q42","q43")
#dimnames(results.array)<- list(1:n.step, parnames, as.character(i.vec))
#dimnames(results.clone)<- list(1:n.step, parnames, as.character(name.vec))


### Writing into the files
write.csv(results.clone,file=paste("clone_taxa", n.tree ,"_",k,"_",mytimeHisse,"_",n.step,".csv",sep = "")) 
# clone result
