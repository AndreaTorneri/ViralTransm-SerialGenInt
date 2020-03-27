#####################################################################
#
# simFunctionsSerGen.R  
#
#
# This script contains a collection of functions used to simulate an epidemic dynamic based on the notion 
# of effective contact process. The Infection dynamic is described by contact rates accepted/rejected according to an infectiousness measure
# based on the viral load. Control strategies that decrease infectiousness or contact rate are included in the script.
# Precisely, we simulate:
# Isolation - infected individuals in isolation are assumed to do not make contacts anymore, i.e. contact rate is set to zero
# Quarantine - infected individuals in isolation are assumed to decrease their contact rate
# Antiviral drug - infected individuals injected with an antiviral decrease their infectiousness measure
#
# In the first part there are functions used in the simulators, as the simulation of the incubation period, infectiousness measure
# and symptomatic period. In a second part we report the simulators used to simulate the epidemic outbreak 
#
# For the function reported hereunder is specified the input and the output parameters.
# 
#
# Author: Andrea Torneri
# Last version: 20/03/2020
#############################################################################################

#############################################################################################
# nCov.diagnosis: Function to simulate the time to diagnosis (Exponentially distributed)
#
# The time to diagnosis is intended as the time from symptoms onset to the diagnosis.
#
#INPUT:
# avg.dt - average time to diagnosis (in day)
#
#OUTPUT:
# t.d - time to diagnosis
#


nCov.diagnosis<-function(avg.dt,distr.dt,sd.dt){
  if (distr.dt=="const"){
    t.d<-avg.dt
  }
  if (distr.dt=="exp"){
    t.d<-rexp(1,rate = 1/avg.dt) 
  }
  if (distr.dt=="gam"){
    sc.<-sd.dt^2/avg.dt
    sh.<-(avg.dt/sd.dt)^2
    t.d<-rgamma(1, shape = sh., scale = sc.)
  }
  return(t.d)
}


#############################################################################################
# nCov.IncubationPeriod: Function to simulate the Incubation Period
#
# The time to incubation is the time to symptom onset since infection
# This is here simulated via a Weibull distribution with shape 2. This because, together with 
# the mean value of 5.2, this represents the incubation period for COVID-19 (Zhang et al. 2020 MedRxiv).
# The function can be easily modified for other distributions
#
#
#INPUT:
# avg.inc - average incubation time (in day)
#
#OUTPUT:
#  - incubation period
#

nCov.IncubationPeriod<-function(avg.inc){
  shape<-2
  scale<-avg.inc*shape/gamma(1/shape)
  return(rweibull(1,shape = shape, scale = scale))
}


######################################################################################################
# nCov.InfMeasure: Function evaluate the infectiousness measure at a precise time point since symptoms onset
#
# The infectiousness measure is set to have a similar shape to the viral load of the specific pathogen.
# Here we want to represent the viral load for COVID-19. Therefore we set the curve to represent the few
# data available to date. (Pan et al. 2020, Lancet; Zhou et al. 2020 Nature; Zou et al. 2020 NEJM)
#
# We choose to use a gamma distribution. When more information would be available, the script can be easily adapted
# The infectiousness measure is scaled to have the same shape for different length of the infectious period.
# The curve is defined over the incubation and symptomatic period, or alternatively, over the exposed and infectious period.
#
#
# INPUT:
# t - time since infection
# lengthI - length of the symptomatic plus incubation period 
#
# OUTPUT:
# value of the infectiousness measure at the specific time since infection
# 

nCov.InfMeasure<-function(t, lengthI){
  vload.comp<-(25/lengthI)*dgamma(25/lengthI*t, shape =7 , scale =1.43 )
  return(vload.comp)
  
}

######################################################################################################
# nCov.SymptmPeriod: Function to simulate the infectious period length
#
# This is modeled to follow an exponential distribution.
#
# INPUT:
# mu.IP: mean length of the symptomatic period
#
# OUTPUT:
# continous value for the length of the symptomatic period
#   

#Function to simulate Infectious period. Assumed to be exponentially distributed
nCov.SymptmPeriod<-function(mu.IP){
  return(rexp(1,1/mu.IP))
  
  #return(rgamma(1,shape = 37.60435, scale = 0.2473118))
  
}




nCov.ctaRate<-function(R0,k.overd,distr.R0,lambda){
  if (distr.R0=="Nbin"){
    nctr<-rnbinom(1,size = k.overd,mu=R0)*lambda/R0
  }
  else{
    nctr<-lambda
  }
  return(nctr)
}

####################################################################################
####################################################################################
#Functions to compute the mean realized generation and serial interval over an epidemic outbreak.
# This quantinty is computed as the average of all the serial interval, or generation interval, realized
# in an entire simulation.
#
# Functions:
#
# nCov.FWDGenerationTime : compute the forward generation time over time.
# nCov.MeanGenerationTime : compute the mean among all the generation times realized.
# nCov.MeanGenerationTime.truncated : compute the mean generation time stopping the observation at a time t

# Similar functions are set for the serial interval.

nCov.FWDGenerationTime<-function(status.matrix){
  out<-matrix(NA,1,2)
  temp.inf<-0
  name<-c("time","mean fwd")
  dimnames(out)<-list(NULL,name)
  inf.time<-sort(status.matrix[!is.na(status.matrix[,2]),2]) #time of infection
  
  #initializing variables
  if (length(inf.time)>1){
    for(i in inf.time){
      infector<-which(status.matrix[,2]==i)
      if(length(infector)>1){ #if there are more than one individual infected at the same time-point. This is the case when there is an initial multiple seeding
        temp.inf<-infector
        temp.gt<-0
        for (n in 1:length(temp.inf)){
          infector<-temp.inf[n]
          infectee<- which(status.matrix[,3]==infector)
          if (length(infectee)>0){ #if the infector infects someone
            temp<-0
            for(j in 1:length(infectee)){
              temp[j]<-status.matrix[infectee[j],2]-status.matrix[infector,2]
            }
          }
          temp.gt[n]<-mean(temp)
        }
        
        out<-rbind(out,c(inf.time[i],mean(temp.gt)))
      }else{
        infectee<- which(status.matrix[,3]==infector)
        if (length(infectee)>0){ #if the infector infects someone
          temp<-0
          for(j in 1:length(infectee)){
            temp[j]<-status.matrix[infectee[j],2]-status.matrix[infector,2]
          }
        }
        out<-rbind(out,c(inf.time[i],mean(temp)))
      }
    }
  }
  
  return(out)
  
}

nCov.MeanGenerationTime<-function(status.matrix){
  Tg<-NULL
  temp.inf<-0
  inf.time<-sort(status.matrix[!is.na(status.matrix[,2]),2]) #time of infection
  
  #initializing variables
  if (length(inf.time)>1){
    for(i in inf.time){
      infector<-which(status.matrix[,2]==i)
      if(length(infector)>1){ #if there are more than one individual infected at the same time-point. This is the case when there is an initial multiple seeding
        temp.inf<-infector
        temp.gt<-0
        for (n in 1:length(temp.inf)){
          infector<-temp.inf[n]
          infectee<- which(status.matrix[,3]==infector)
          if (length(infectee)>0){ #if the infector infects someone
            for(j in 1:length(infectee)){
              Tg<-c(Tg,status.matrix[infectee[j],2]-status.matrix[infector,2])
            }
          }
        }
      }else{
        infectee<- which(status.matrix[,3]==infector)
        if (length(infectee)>0){ #if the infector infects someone
          for(j in 1:length(infectee)){
            Tg<-c(Tg,status.matrix[infectee[j],2]-status.matrix[infector,2])
          }
        }
      }
    }
  }
  
  ifelse(is.null(Tg),Tg<-NA,Tg<-Tg)
  
  return(mean(Tg))
  
}

nCov.MeanGenerationTime.truncated<-function(status.matrix,t){
  Tg<-NULL
  temp.inf<-0
  inf.time<-sort(status.matrix[!is.na(status.matrix[,2]),2]) #time of infection
  inf.time<-inf.time[which(inf.time<t)]
  
  #initializing variables
  if (length(inf.time)>1){
    for(i in inf.time){
      infector<-which(status.matrix[,2]==i)
      if(length(infector)>1){ #if there are more than one individual infected at the same time-point. This is the case when there is an initial multiple seeding
        temp.inf<-infector
        temp.gt<-0
        for (n in 1:length(temp.inf)){
          infector<-temp.inf[n]
          infectee<- which(status.matrix[,3]==infector)
          if (length(infectee)>0){ #if the infector infects someone
            for(j in 1:length(infectee)){
              if(status.matrix[infectee[j],2]<t){
                Tg<-c(Tg,status.matrix[infectee[j],2]-status.matrix[infector,2])
              }
              
            }
          }
        }
      }else{
        infectee<- which(status.matrix[,3]==infector)
        if (length(infectee)>0){ #if the infector infects someone
          for(j in 1:length(infectee)){
            if(status.matrix[infectee[j],2]<t){
              Tg<-c(Tg,status.matrix[infectee[j],2]-status.matrix[infector,2])
            }
          }
        }
      }
    }
  }
  
  ifelse(is.null(Tg),Tg<-NA,Tg<-Tg)
  
  return(mean(Tg))
  
}

nCov.FWDSerialInterval<-function(status.matrix){
  out<-matrix(NA,1,2)
  temp.inf<-0
  name<-c("time","mean fwd")
  dimnames(out)<-list(NULL,name)
  symptomatic<-c(which(status.matrix[,4]==1),which(status.matrix[,4]==2))
  inf.time<-sort(status.matrix[symptomatic,2]) #time of infection
  
  #initializing variables
  if (length(inf.time)>1){
    for(i in inf.time){
      infector<-which(status.matrix[,2]==i)
      if(length(infector)>1){ #if there are more than one individual infected at the same time-point. This is the case when there is an initial multiple seeding
        temp.inf<-infector
        temp.gt<-NULL
        for (n in 1:length(temp.inf)){
          infector<-temp.inf[n]
          infectee<- which(status.matrix[,3]==infector)
          if (length(infectee)>0){ #if the infector infects someone
            temp<-NULL
            for(j in 1:length(infectee)){
              if (status.matrix[infectee[j],4]!=0){
                temp[j]<-status.matrix[infectee[j],5]-status.matrix[infector,5]
              }
            }
          }
          if (!is.null(temp)){
            temp.gt<-c(temp.gt, mean(temp)) 
          }
        }
        if (!is.null(temp.gt)){out<-rbind(out,c(i,mean(temp.gt)))}
      }else{
        infectee<- which(status.matrix[,3]==infector)
        if (length(infectee)>0){ #if the infector infects someone
          temp<-NULL
          for(j in 1:length(infectee)){
            if (status.matrix[infectee[j],4]!=0){
              temp[j]<-status.matrix[infectee[j],5]-status.matrix[infector,5]
            }
          }
        }
        if (!is.null(temp)){out<-rbind(out,c(i,mean(temp)))}
      }
    }
  }
  
  return(out)
  
}

nCov.MeanSerialInterval<-function(status.matrix){
  Tg<-NULL
  temp.inf<-0
  symptomatic<-c(which(status.matrix[,4]==1),which(status.matrix[,4]==2))
  inf.time<-sort(status.matrix[symptomatic,2]) #time of infection
  
  
  #initializing variables
  if (length(inf.time)>1){
    for(i in inf.time){
      infector<-which(status.matrix[,2]==i)
      if(length(infector)>1){ #if there are more than one individual infected at the same time-point. This is the case when there is an initial multiple seeding
        temp.inf<-infector
        temp.gt<-0
        for (n in 1:length(temp.inf)){
          infector<-temp.inf[n]
          infectee<- which(status.matrix[,3]==infector)
          if (length(infectee)>0){ #if the infector infects someone
            for(j in 1:length(infectee)){
              if (status.matrix[infectee,4]!=0){
                Tg<-c(Tg,status.matrix[infectee[j],5]-status.matrix[infector,5])
              }
            }
          }
        }
      }else{
        infectee<- which(status.matrix[,3]==infector)
        if (length(infectee)>0){ #if the infector infects someone
          for(j in 1:length(infectee)){
            if (status.matrix[infectee[j],4]!=0){
              Tg<-c(Tg,status.matrix[infectee[j],5]-status.matrix[infector,5])
            }
          }
        }
      }
    }
  }
  
  ifelse(is.null(Tg),mean.Tg<-NA,mean.Tg<-mean(Tg))
  
  return(mean.Tg)
  
  
}

nCov.MeanSerialInterval.truncated<-function(status.matrix,t){
  Tg<-NULL
  temp.inf<-0
  symptomatic<-c(which(status.matrix[,4]==1),which(status.matrix[,4]==2))
  inf.time<-sort(status.matrix[symptomatic,2]) #time of infection
  inf.time<-inf.time[which(inf.time<t)]
  
  
  #initializing variables
  if (length(inf.time)>1){
    for(i in inf.time){
      infector<-which(status.matrix[,2]==i)
      if(length(infector)>1){ #if there are more than one individual infected at the same time-point. This is the case when there is an initial multiple seeding
        temp.inf<-infector
        temp.gt<-0
        for (n in 1:length(temp.inf)){
          infector<-temp.inf[n]
          infectee<- which(status.matrix[,3]==infector)
          if (length(infectee)>0){ #if the infector infects someone
            for(j in 1:length(infectee)){
              if (status.matrix[infectee[j],4]!=0 & status.matrix[infectee[j],2]<t ){
                Tg<-c(Tg,status.matrix[infectee[j],5]-status.matrix[infector,5])
              }
            }
          }
        }
      }else{
        infectee<- which(status.matrix[,3]==infector)
        if (length(infectee)>0){ #if the infector infects someone
          for(j in 1:length(infectee)){
            if (status.matrix[infectee[j],4]!=0 & status.matrix[infectee[j],2]<t){
              Tg<-c(Tg,status.matrix[infectee[j],5]-status.matrix[infector,5])
            }
          }
        }
      }
    }
  }
  
  ifelse(is.null(Tg),mean.Tg<-NA,mean.Tg<-mean(Tg))
  
  return(mean.Tg)
  
  
}



####################################################################################
####################################################################################
# Simulators
# Here we report the scripts used to simulate the epidemic outbreak. Scripts are different based on the control strategy 
# they represent. More precisely:
#
#
#
# We report here the list of input parameters needed in the simulator
# n : population size
# lambda: contact rate
# sigma: probability of severe symptoms
# avg.inc: mean incubation period 
# avg.dt: mean time to diagnosis
# eta: probability of succesful tracing 
# rate.quar: contact rate in quarantine
# 



# OUTPUT
# Time events: matrix with three columns indicating the time at which an infection/recover event happens (first column) that time at which this event happens (second column) and the individual id (third column)
# Status.matrix: matrix that keep track of different features of the individaul during the epidemic, e.g. status (infected/recovered/susceptible), time of infection, time of symptom onsets, time of diagnosis, time of treatment.




nCov.simulator.QuarIso<-function(n, lambda, rho, sigma,avg.inc,R.a,R.s,R.ss,rq,ri, avg.dt, distr.dt, sd.dt,distr.R0,k.overd, mu.a,mu.s,mu.ss ){
  
  status.matrix <- matrix(NA,nrow = n,ncol = 6) #matrix containing information about the state of the individuals
  col.name<-c("infected","time.of.infection","infector", "severity", "TimeSymptomOnset","TimeQuarantine")
  dimnames(status.matrix)<-list(NULL,col.name)
  status.matrix[,1]<-0 #all the population is initially susceptible
  recovery.vector<-rep(NA,n) #vector of recovery times
  quarantine.day<-rep(Inf,n) #day at which individuals will be quarantined
  infectives<-rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
  quarantine<-rep(0,n) #vector that says who is at in quarantine at the current time. 0 means no quarantine, 1 quarantine
  current.time<-0
  index.contact<-rep(0,n) # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
  transmission.parameters<-data.frame("id"=1:n,"q"=rep(NA,n),"total_infectionPeriod"=rep(NA,n), "contact_rate"=lambda)   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  # first infected: randomly chosen in the population (among the susceptibles)
  first<-sample(which(status.matrix[,1]==0), 1) #initial case
  status.matrix[first,1] <- 1 
  status.matrix[first,2] <- 0
  status.matrix[first,5]<-current.time+nCov.IncubationPeriod(avg.inc =avg.inc) #(Exposed period - this is assume to be of same length of the incubation period)
  infectives[first]<-1
  time.events<-matrix(NA,1,3)
  if (runif(1)<rho){ #check if the individual shows symtpoms
    if (runif(1)<sigma){ #check whether are severe
      status.matrix[first,4]<-2
      recovery.vector[first]<-status.matrix[first,5]+nCov.SymptmPeriod(mu.IP = mu.ss) # the total length since infection (Exposed+IP) 
      transmission.parameters$q[first]<-R.ss/lambda
      transmission.parameters$contact_rate[first]<-nCov.ctaRate(R0=R.ss, distr.R0 = distr.R0, k.overd = k.overd, lambda = lambda)
      if (transmission.parameters$contact_rate[first]>0){
        index.contact[first]<-1 #when 1 individual proposes a contact 
      }
      quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt, distr.dt = distr.dt, sd.dt = sd.dt)
      time.events[1,]<-c(current.time,1.2,first)
      status.matrix[first,6]<-quarantine.day[first]
    }else{
      status.matrix[first,4]<-1
      recovery.vector[first]<-status.matrix[first,5]+nCov.SymptmPeriod(mu.IP = mu.s) # the total length since infection (Exposed+IP) 
      transmission.parameters$q[first]<-R.s/lambda
      transmission.parameters$contact_rate[first]<-nCov.ctaRate(R0=R.s, distr.R0 = distr.R0, k.overd = k.overd, lambda = lambda)
      if (transmission.parameters$contact_rate[first]>0){
        index.contact[first]<-1 #when 1 individual proposes a contact 
      }
      quarantine.day[first]<-status.matrix[first,5]+nCov.diagnosis(avg.dt = avg.dt, distr.dt = distr.dt, sd.dt = sd.dt)
      status.matrix[first,6]<-quarantine.day[first]
      time.events[1,]<-c(current.time,1.1,first)
    }
  }else{
    status.matrix[first,4]<-0
    recovery.vector[first]<-status.matrix[first,5]+nCov.SymptmPeriod(mu.IP = mu.a) # the total length since infection (Exposed+IP) 
    transmission.parameters$q[first]<-R.a/lambda
    transmission.parameters$contact_rate[first]<-nCov.ctaRate(R0=R.a, distr.R0 = distr.R0, k.overd = k.overd, lambda = lambda)
    if (transmission.parameters$contact_rate[first]>0){
      index.contact[first]<-1 #when 1 individual proposes a contact 
    }
    time.events[1,]<-c(current.time,1.0,first)
  }
  transmission.parameters$total_infectionPeriod[first]<-recovery.vector[first]-current.time
  contact.time<-data.frame("id"=1:n,"pr.ctc"=rep(NA,n))   #matrix containing the proposed time of the next contact (first colum) and the contact individual (second column)
  proposed.individual<-0
  int.time<-0
  temp.contact.time<-0
  indiv.prop.ctc<-0
  T_NextCtc<-0
  recovered<-0
  
  #When only the first pathogen is present
  while(sum(infectives, na.rm = TRUE) > 0){ #while there are still infectives
    #Phase 1: individuals who has to, propose a new social contact
    for (i in which(index.contact==1) ){ # for all the individuals that has to propose a global contact
      temp.contact.time<-rexp(1,transmission.parameters$contact_rate[i])+current.time# I generate the next interarrival time for individual i
      index.contact[i]<-0
      contact.time$pr.ctc[i]<-temp.contact.time # If it is an infectious contact I keep track of the individual and of the time
    }
    
    #Phase 2: identify the next event: possible infection, someone is quarantined or a recovery
    ifelse(length(which(is.na(contact.time$pr.ctc)==FALSE))>0,T_NextCtc<-min(contact.time$pr.ctc, na.rm = T),T_NextCtc<-Inf) # among all the proposed social contact between houeholds we select the minimum
    ifelse(length(which(!is.infinite(quarantine.day)))>0,Q<-min(quarantine.day),Q<-Inf ) #minimum quarantine pathogen 1
    R_a<-min(recovery.vector, na.rm = T) # minimum among the recovery times
    
    
    if (Q<=min(T_NextCtc,R_a)){
      current.time<-Q
      person.quarantined<-which(Q==quarantine.day)
      quarantine[person.quarantined]<-1
      quarantine.day[person.quarantined]<-Inf
      if (status.matrix[person.quarantined,4]==2){
       # contact.time$pr.ctc[person.quarantined]<-NA # !!!!!!!!!!
      if (ri==0){
        contact.time$pr.ctc[person.quarantined]<-NA
      }else{
       transmission.parameters$contact_rate[person.quarantined]<-ri
      }
      }else{
        transmission.parameters$contact_rate[person.quarantined]<-rq
      }
    }else{
      if (T_NextCtc<R_a){
        current.time<-T_NextCtc
        infector<-which(contact.time$pr.ctc ==T_NextCtc) 
        infectee<-sample(setdiff(1:n,infector),1)
        acc.rate<-nCov.InfMeasure(t=(current.time-status.matrix[infector,2]), lengthI = transmission.parameters$total_infectionPeriod[infector])*transmission.parameters$q[infector]
        if (quarantine[infectee]==1){acc.rate<-0}
        if (status.matrix[infectee,1]==0 & runif(1)<acc.rate){
          status.matrix[infectee,1]<-1
          status.matrix[infectee,2]<-current.time
          status.matrix[infectee,3]<-infector
          status.matrix[infectee,5]<-current.time+nCov.IncubationPeriod(avg.inc=avg.inc) #(Exposed period - this is assume to be of same length of the incubation period)
          infectives[infectee]<-1
          if (runif(1)<rho){ #check if the individual shows symtpoms
            if (runif(1)<sigma){ #check whether are severe
              status.matrix[infectee,4]<-2
              recovery.vector[infectee]<-status.matrix[infectee,5]+nCov.SymptmPeriod(mu.IP = mu.ss) # the total length since infection (Exposed+IP) 
              transmission.parameters$q[infectee]<-R.ss/lambda
              transmission.parameters$contact_rate[infectee]<-nCov.ctaRate(R0=R.ss, distr.R0 = distr.R0, k.overd = k.overd, lambda = lambda)
              if (transmission.parameters$contact_rate[infectee]>0){
                index.contact[infectee]<-1 #when 1 individual proposes a contact 
              }
              quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt, distr.dt = distr.dt, sd.dt = sd.dt) 
              status.matrix[infectee,6]<-quarantine.day[infectee]
              time.events<-rbind(time.events,c(current.time,1.2,infectee))
            }else{
              status.matrix[infectee,4]<-1
              recovery.vector[infectee]<-status.matrix[infectee,5]+nCov.SymptmPeriod(mu.IP = mu.s) # the total length since infection (Exposed+IP) 
              transmission.parameters$q[infectee]<-R.s/lambda
              transmission.parameters$contact_rate[infectee]<-nCov.ctaRate(R0=R.s, distr.R0 = distr.R0, k.overd = k.overd, lambda = lambda)
              quarantine.day[infectee]<-status.matrix[infectee,5]+nCov.diagnosis(avg.dt = avg.dt, distr.dt = distr.dt, sd.dt = sd.dt) # quarantined after some day from symptom onset              
              status.matrix[infectee,6]<-quarantine.day[infectee]
              if (transmission.parameters$contact_rate[infectee]>0){
                index.contact[infectee]<-1 #when 1 individual proposes a contact 
              }
              time.events<-rbind(time.events,c(current.time,1.1,infectee))
            }
          }else{
            status.matrix[infectee,4]<-0
            recovery.vector[infectee]<-status.matrix[infectee,5]+nCov.SymptmPeriod(mu.IP = mu.a) # the total length since infection (Exposed+IP) 
            transmission.parameters$q[infectee]<-R.a/lambda
            transmission.parameters$contact_rate[infectee]<-nCov.ctaRate(R0=R.a, distr.R0 = distr.R0, k.overd = k.overd, lambda = lambda)
              if (transmission.parameters$contact_rate[infectee]>0){
                index.contact[infectee]<-1 #when 1 individual proposes a contact 
              }
          }
          transmission.parameters$total_infectionPeriod[infectee]<-recovery.vector[infectee]-current.time
          }
        index.contact[infector]<-1
        contact.time$pr.ctc[infector]<-NA
        #Phase 2.3 a recovery occurs
      }else{
        current.time<-R_a
        recovered<-which(recovery.vector==R_a)
        recovery.vector[recovered]<-NA
        status.matrix[recovered,1]<--1
        time.events<-rbind(time.events,c(current.time,0.1,recovered))
        contact.time[recovered,2:3]<-rep(NA,4)
        infectives[recovered]<-NA
        quarantine.day[recovered]<-Inf
        quarantine[recovered]<-0
      }
    }
  }
  #When also the other pathogen is present.
  timev.name<-c("time","event","who")
  dimnames(time.events)<-list(NULL,timev.name)
  
  return(list(time.events=time.events, status.matrix=status.matrix))
}



