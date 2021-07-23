rm(list=ls())
all.lib<-c("ape","diversitree","phytools", "hisse","ggplot2","gridExtra","ggpubr","numDeriv","reshape2","ggthemes","scales") # add libraries here
lapply(all.lib,require,character.only=TRUE)

#colorpalette-black
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#CC6666", "#9999CC", "#66CC99","#999999")

#Cetacean dataset only
#set the working directory and parameters

setwd("/Users/nhans/DiversificationPaper/Diversification2020")
n.step<-15000
k.vec<- 2^seq(1,9,by=1)
sizeofk.vec<-length(k.vec)

# Thinning the MCMC chain
thin<-sequence((n.step-3000)/10,1,10) # thinning interval - 10
l.thin<-length(thin)

#parameters for each model
parnames<-c("tau1","tau2", "tau3", "tau4", "ep1","ep2", "ep3", "ep4", "q12","q13","q14","q21","q23","q24","q31","q32","q34","q41","q42","q43")

 
 
 # HiSSE plots together (1)

##No clones
his.no.clone<-read.csv("/Users/nhans/DiversificationPaper/Diversification2020/posteriorsMCMCwhales_Jul_15_21_50000.csv")[3001:15000,3:22]
oldparnames<-colnames(his.no.clone) # Because i forgot to re-parameterize in the code but named it --- so fix the code on cluster to avoid having to do
colnames(his.no.clone)<-parnames

# Getting the dimension
re.row<-dim(his.no.clone)[1]

# Thinning the MCMC chain
thin<-sequence(re.row/10,1,10) # thinning interval - 10
l.thin<-length(thin)

thin.his.no.clone<-matrix(nrow = l.thin, ncol = 8)
# thinning out the posterior (no clone)
count<-0
for (t in thin){
 count<-count+1
 thin.his.no.clone[count,1]<-eval(parse(text=paste("his.no.clone","$tau1",sep="")))[t]
 thin.his.no.clone[count,2]<-eval(parse(text=paste("his.no.clone","$tau2",sep="")))[t]
 thin.his.no.clone[count,3]<-eval(parse(text=paste("his.no.clone","$tau3",sep="")))[t]
 thin.his.no.clone[count,4]<-eval(parse(text=paste("his.no.clone","$tau4",sep="")))[t]
 
 thin.his.no.clone[count,5]<-eval(parse(text=paste("his.no.clone","$ep1",sep="")))[t]
 thin.his.no.clone[count,6]<-eval(parse(text=paste("his.no.clone","$ep2",sep="")))[t]
 thin.his.no.clone[count,7]<-eval(parse(text=paste("his.no.clone","$ep3",sep="")))[t]
 thin.his.no.clone[count,8]<-eval(parse(text=paste("his.no.clone","$ep4",sep="")))[t]
 
 #thin.his.no.clone[count,9]<-eval(parse(text=paste("his.no.clone","$q12",sep="")))[t]
 #thin.his.no.clone[count,10]<-eval(parse(text=paste("his.no.clone","$q13",sep="")))[t]
 #thin.his.no.clone[count,11]<-eval(parse(text=paste("his.no.clone","$q14",sep="")))[t]
 #thin.his.no.clone[count,12]<-eval(parse(text=paste("his.no.clone","$q21",sep="")))[t]
 
 #thin.his.no.clone[count,13]<-eval(parse(text=paste("his.no.clone","$q23",sep="")))[t]
 #thin.his.no.clone[count,14]<-eval(parse(text=paste("his.no.clone","$q24",sep="")))[t]
 #thin.his.no.clone[count,15]<-eval(parse(text=paste("his.no.clone","$q31",sep="")))[t]
 #thin.his.no.clone[count,16]<-eval(parse(text=paste("his.no.clone","$q32",sep="")))[t]
 
 #thin.his.no.clone[count,17]<-eval(parse(text=paste("his.no.clone","$q34",sep="")))[t]
 #thin.his.no.clone[count,18]<-eval(parse(text=paste("his.no.clone","$q41",sep="")))[t]
 #thin.his.no.clone[count,19]<-eval(parse(text=paste("his.no.clone","$q42",sep="")))[t]
 #thin.his.no.clone[count,20]<-eval(parse(text=paste("his.no.clone","$q43",sep="")))[t]
 
}

eigen.max.his<-max(eigen(cov(thin.his.no.clone))$values)
eigen.max.his
# 30.83422 with transition rates
# 0.04000248 without transition rates 

# Clones 
# Load and discard 20% burnin
for (k in k.vec){
  temp.dat<-read.csv(paste("/Users/nhans/DiversificationPaper/Diversification2020/clone_taxa_MCMCwhales_",k,"_Jul_15_21_50000.csv",sep=""))[3001:15000,2:21]
  colnames(temp.dat)<-parnames
  assign(paste("his.clone.cet",k,sep=""), temp.dat)
}

# Resetting dimensions
re.row<-dim(his.clone.cet2)[1]

# Thinning the MCMC chain
thin<-sequence(re.row/10,1,10) # thinning interval - 10
l.thin<-length(thin)

#Format: 1=0A, 2=1A, 3=0B, 4=1B.
# Assign empty matrix for 4 turnover rates and 4 extinction fractions 
for (i in 1:4){
  temp.mat<-matrix(NA, ncol=sizeofk.vec, nrow=l.thin,dimnames=list(as.character(1:l.thin),as.character(k.vec)))
  assign(paste("tau",i, ".his.clone",sep = ""), temp.mat) # turnover rates
  assign(paste("ep",i, ".his.clone",sep = ""), temp.mat)  # extinction fraction
  assign(paste("trans",i, ".his.clone",sep = ""), temp.mat)  # transition rates (note that there are 12)
}


#------------------------------------------Saving turnover rates, extinction fractions in different matrix
counter<-0
for(k in k.vec)
{
  counter<-counter+1
  rowcount<-0
  for(t in thin){
  rowcount<-rowcount+1
  tau1.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$tau1",sep="")))[t]
  tau2.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$tau2",sep="")))[t]
  tau3.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$tau3",sep="")))[t]
  tau4.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$tau4",sep="")))[t]

  ep1.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$ep1",sep="")))[t]
  ep2.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$ep2",sep="")))[t]
  ep3.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$ep3",sep="")))[t]
  ep4.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$ep4",sep="")))[t]
  
  # Transition rates
  trans1.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$q12",sep="")))[t]
  trans2.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$q13",sep="")))[t]
  trans3.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$q14",sep="")))[t]
  trans4.his.clone[rowcount,counter]<-eval(parse(text=paste("his.clone.cet",k,"$q21",sep="")))[t]
 }
}


#----------------------------------------------- Reparameterizing from tau and ep to lambda and mu
# lambda1<-tau1/(1+ep1)
# mu1<- lambda1* ep1
#Reparameterizing
lambda1.his.clone<-tau1.his.clone/(1+ep1.his.clone)
lambda2.his.clone<-tau2.his.clone/(1+ep2.his.clone)
lambda3.his.clone<-tau3.his.clone/(1+ep3.his.clone)
lambda4.his.clone<-tau4.his.clone/(1+ep4.his.clone)

mu1.his.clone<- lambda1.his.clone* ep1.his.clone
mu2.his.clone<- lambda2.his.clone* ep2.his.clone
mu3.his.clone<- lambda3.his.clone* ep3.his.clone
mu4.his.clone<- lambda4.his.clone* ep4.his.clone




