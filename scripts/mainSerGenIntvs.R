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
cat(",R.s=",R.s)
mu.a = as.numeric(args[12]) # reproduction number asymptomatic
cat(",mu.a=",mu.a)
mu.s = as.numeric(args[13]) # reproduction number symptomatic
cat(",mu.s=",mu.s)
mu.ss = as.numeric(args[14]) # reproduction number symptomatic
cat(",mu.ss=",mu.ss)
scenario = args[15] # reproduction number symptomatic
cat(",scenario=",scenario)
k.overd = as.numeric(args[16]) # reproduction number symptomatic
cat(",mu.ss=",k.overd)


#load packages
library(foreach)
library(doParallel)
# load functions
source("SimFunctionsSerGen.R")
# Parameters
n<-100 #Population size
lambda<-12 #daily rate of social contacts
sigma<- 0.16 # proportion of severe cases 
avg.inc<-5.2 # NEJM Li
distr.R0<-"Nbin"
lambda.q<-lambda*rq
lambda.i<-lambda*ri

cl<-makeCluster(as.numeric(cores))
registerDoParallel(cl)

nSim<-50
epi.outbreak<-list()
nSeed<-14022020
set.seed(nSeed)

epi.outbreak<-foreach(i = 1:nSim) %dopar%{
    epi.outbreak[[i]]<-nCov.simulator.QuarIso(n=n,lambda = lambda, rho = rho, sigma = sigma, avg.inc = avg.inc, R.a = R.a, R.s = R.s, R.ss=R.ss, rq=lambda.q, ri=lambda.i, avg.dt = avg.dt, distr.dt = distr.dt, sd.dt = sd.dt, distr.R0 = distr.R0,k.overd = k.overd,mu.a = mu.a, mu.s = mu.s, mu.ss = mu.ss)
}

stopCluster(cl)

finalSize<-NULL
not.extinct<-NULL
for (i in 1:nSim){
  finalSize[i]<-length(c(which(epi.outbreak[[i]]$time.events[,2]==1.0),which(epi.outbreak[[i]]$time.events[,2]==1.1),which(epi.outbreak[[i]]$time.events[,2]==1.2)))
  if (finalSize[i]>round(n*0.1)){not.extinct<-c(not.extinct,i)}
}

MeanSerialInterval<-NULL
MeanGenTime<-NULL

for (i in not.extinct) {
  MeanSerialInterval<-c(MeanSerialInterval, nCov.MeanSerialInterval(status.matrix=epi.outbreak[[i]]$status.matrix)) 
  MeanGenTime<-c(MeanGenTime, nCov.MeanGenerationTime(status.matrix=epi.outbreak[[i]]$status.matrix)) 
}

index.of.generations<-NULL
infectee<-0
infector<-0
prop.before<-NULL
prop.after<-NULL
for (i in not.extinct){
  gen.before.symptoms<-0
  gen.after<-0
  timEv<-epi.outbreak[[i]]$time.events
  StMat<-epi.outbreak[[i]]$status.matrix
  index.of.generations<-sort(c(which(timEv[,2]==1.2),which(timEv[,2]==1.1)))
  temp.gen<-NULL
  for (j in 2:length(index.of.generations)){
    infectee<-timEv[index.of.generations[j],3]
    infector<-StMat[infectee,3]
    if (StMat[infectee,2]<StMat[infector,5]){
      gen.before.symptoms<-gen.before.symptoms+1
    }else{
      gen.after<-gen.after+1
    }
  }
  prop.before<-c(prop.before,gen.before.symptoms/(gen.before.symptoms+gen.after))
  prop.after<-c(prop.after,gen.after/(gen.before.symptoms+gen.after))
}


name<-paste("SimNCovidGenSER",n,"_Rs",R.s,"_Ra",R.a,"_Scenario",scenario,"_lq",lambda.q,"_td",avg.dt,"_tdDistr",distr.dt,"_muA",mu.a,"_muS",mu.s,"_.RData", sep = "")

save(epi.outbreak, MeanGenTime, MeanSerialInterval, prop.before, file = name)

