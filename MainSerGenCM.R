args <- commandArgs(trailingOnly = TRUE)

out = args[1] #working directory
cat(",out=",out)
cores = as.numeric(args[2]) #number of cores to run in parallel
cat(",cores=",cores)
avg.dt = as.numeric(args[3]) #number of cores to run in parallel
cat(",avg.dt=",avg.dt)
distr.dt = args[4] #scenario: const - exp - gam
cat(",distr.dt=",distr.dt)
sd.dt = as.numeric(args[5]) # st deviation distr.dt (in case gamma distributed, otherwise set 0)
cat(",sd.dt=",sd.dt)
rho = as.numeric(args[6]) #proportion of diagnosed
cat(",rho=",rho)
rq = as.numeric(args[7]) # decrease in contact rate due to quarantine
cat(",rq=",rq)
ri = as.numeric(args[8]) # decrease in contact rate due to quarantine
cat(",ri=",ri)
R.a = as.numeric(args[9]) # reproduction number asymptomatic
cat(",R.a=",R.a)
R.s = as.numeric(args[10]) # reproduction number symptomatic
cat(",R.s=",R.s)
R.ss = as.numeric(args[11]) # reproduction number symptomatic
cat(",R.s=",R.ss)
mu.a = as.numeric(args[12]) # reproduction number asymptomatic
cat(",mu.a=",mu.a)
mu.s = as.numeric(args[13]) # reproduction number symptomatic
cat(",mu.s=",mu.s)
mu.ss = as.numeric(args[14]) # reproduction number symptomatic
cat(",mu.ss=",mu.ss)
scenario = args[15] # reproduction number symptomatic
cat(",scenario=",scenario)
k.overd = as.numeric(args[16]) # reproduction number symptomatic
cat(",k.overd=",k.overd)
distr.R0 = args[17] # reproduction number symptomatic
cat(",distr.R0=",distr.R0)
InfM = args[18] # reproduction number symptomatic
cat(",InfM=",InfM)
day.sympt.peak=args[19]
cat(",day.sympt.peak=",day.sympt.peak)




print(paste0(out,cores,avg.dt,distr.dt,sd.dt, rho, rq, ri, R.a, R.s, R.ss, mu.a, mu.s, mu.ss,scenario,k.overd,distr.R0,InfM,day.sympt.peak))


#load packages
library(foreach)
library(doParallel)
library(doRNG)
#library(sn)
# load functions
source("SimFunctionsSerGen.R")
# Parameters
n<-100 #Population size
lambda<-12 #daily rate of social contacts
sigma<- 0.16 # proportion of severe cases 
avg.inc<-5.2 # NEJM Li
lambda.q<-lambda*rq
lambda.i<-lambda*ri

cl<-makeCluster(as.numeric(cores))
registerDoParallel(cl)

nSim<-100000
epi.outbreak<-list()
nSeed<-14022020
set.seed(nSeed)
symptInf.val<-0.12

  epi.outbreak<-foreach(i = 1:nSim) %dopar%{
    epi.outbreak[[i]]<-nCov.simulator.QuarIso.varyingSOV2(n=n,lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, R.a = R.a, R.s = R.s, R.ss=R.ss, rq=lambda.q, ri=lambda.i, avg.dt = avg.dt, distr.dt = distr.dt, sd.dt = sd.dt, distr.R0 = distr.R0,k.overd = k.overd,mu.a = mu.a, mu.s = mu.s, mu.ss = mu.ss, day.sympt.peak = day.sympt.peak)  
  }

stopCluster(cl)

finalSize<-NULL
not.extinct<-NULL
for (i in 1:nSim){
  finalSize[i]<-length(c(which(epi.outbreak[[i]]$time.events[,2]==1.0),which(epi.outbreak[[i]]$time.events[,2]==1.1),which(epi.outbreak[[i]]$time.events[,2]==1.2)))
  if (finalSize[i]>round(n*0.1)){not.extinct<-c(not.extinct,i)}
}
FinSize<-finalSize[not.extinct]

casesPeak<-NULL
for (j in not.extinct) {
  time.events<-epi.outbreak[[j]]$time.events  
  epi.curve<-0
  for (i in 1:length(time.events[,1]) ){
    epi.curve[i]<-length(c(which(time.events[1:i,2]==1.1),which(time.events[1:i,2]==1.2),which(time.events[1:i,2]==1.0)))-length(which(time.events[1:i,2]==0.1))
  }
  casesPeak<-c(casesPeak, max(epi.curve))
}
PeakInc<-casesPeak



MeanSerialInterval<-NULL
SDSerialInterval<-NULL
MeanGenTime<-NULL
SDGenTime<-NULL

for (i in not.extinct) {
  MeanSerialInterval<-c(MeanSerialInterval, nCov.MeanSerialInterval(status.matrix=epi.outbreak[[i]]$status.matrix)$Mean) 
  SDSerialInterval<-c(SDSerialInterval, nCov.MeanSerialInterval(status.matrix=epi.outbreak[[i]]$status.matrix)$Sd) 
  MeanGenTime<-c(MeanGenTime, nCov.MeanGenerationTime(status.matrix=epi.outbreak[[i]]$status.matrix)$Mean) 
  SDGenTime<-c(SDGenTime, nCov.MeanGenerationTime(status.matrix=epi.outbreak[[i]]$status.matrix)$Sd)
}



name<-paste(scenario,n,"_Rs",R.s,"_Ra",R.a,"DistrR",distr.R0,"_muS",mu.s,"_.RData", sep = "")
setwd(out)
save(epi.outbreak,FinSize,PeakInc, MeanSerialInterval,MeanGenTime,SDSerialInterval,SDGenTime, file = name)

