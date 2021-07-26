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
samples.noclone<-mcmc(lik.func, temp.p, nsteps=n.step, prior=temp.prior,lower=0, w=rep(1,20), print.every=0, save.every=100,save.file=paste("posteriorsMCMC", n.tree ,"_",mytimeHisse,"_",n.step,".csv",sep = ""))

count<-0
#Simulations for Posteriors-MCMC
message.NH5<- paste("Saving results:", count , sep= "")
print(message.NH5)
write.csv(samples.noclone,file=paste("posteriors", n.tree ,"_",mytimeHisse,"_",n.step,".csv",sep = "")) 

