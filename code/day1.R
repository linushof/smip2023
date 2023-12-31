pacman::p_load("rtdists", "tidyverse", "msm")


# Sessio 1

# noisy evidence accumulation 

iter <- 1000
for(i in 1:iter){
  # evidence
  e <- sample(c(2, -1), 1, replace = TRUE) + rnorm(1, 0,1)
  v[i + 1] <- v[i] + e
}


# solution

# step 1

simDiffusion=function(v,s,dt,maxiter) {
  x = rep(NA,maxiter+1)
  x[1] = z
  for (i in 1:maxiter) {
    x[i+1] = x[i] + v*dt + s*sqrt(dt)*rnorm(1)
  }
  return(x)
}

v = 1
s = 1
z = 0
dt = 0.01
maxiter = 1000

tmp=simDiffusion(v=v,s=s,dt=dt,maxiter=maxiter)

plot(dt*1:(maxiter+1),tmp,xlab="Time",ylab="Evidence",type="l")


# step 2

simDiffusion=function(v,a,ter,z,s,dt,maxiter) {
  resp = -1
  rt = -1
  x = rep(NA,maxiter+1)
  x[1] = z # starting point
  for (i in 1:maxiter) {
    x[i+1] = x[i] + v*dt + s*sqrt(dt)*rnorm(1)
    if (x[i+1] > a) { # upper threshold
      resp = 2 # which response is made
      rt = ter + i*dt - dt/2 # response time 
      break
    } 
    if (x[i+1] < 0) { # lower threshold
      resp = 1 
      rt = ter + i*dt - dt/2 # response time
      break
    }
  }
  return(list(resp=resp,rt=rt,x=x))
}

v=1
a=2
ter=0.3
z=a/2
s=1
dt=0.01
maxiter=1000

tmp=simDiffusion(v=v,a=a,ter=ter,z=z,s=s,dt=dt,maxiter=maxiter)
tmp

plot(dt*1:(maxiter+1),tmp$x,xlab="Time",ylab="Evidence",type="l",ylim=c(-a/2,a+a/2))
abline(h=a,col="red")
abline(h=0,col="red")

cat("Response =",ifelse(tmp$resp==2,"Correct","Error"),"\n")
cat("Response Time =",tmp$rt,"\n")


# step 3 and 5

simDiffusion=function(v,a,ter,z,s,dt,maxiter,ntrials) {
  
  tdata <- vector("list", length = ntrials)
  
  for(j in 1:ntrials){ 
    
  resp = -1
  rt = -1
  x = rep(NA,maxiter+1)
  x[1] = z # starting point
  for (i in 1:maxiter) {
    x[i+1] = x[i] + v*dt + s*sqrt(dt)*rnorm(1)
    if (x[i+1] > a) { # upper threshold
      resp = 2 # which response is made
      rt = ter + i*dt - dt/2 # response time 
      break
    } 
    if (x[i+1] < 0) { # lower threshold
      resp = 1 
      rt = ter + i*dt - dt/2 # response time
      break
    }
  }
  tdata[[j]] <- list(resp=resp,rt=rt,x=x)
  }
}

v=1
a=2
ter=0.3
z=a/2
s=1
dt=0.01
maxiter=1000
ntrials = 100
simDiffusion(v=v,a=a,ter=ter,z=z,s=s,dt=dt,maxiter=maxiter, ntrials = ntrials)
tdata
tdata[[1]]


# solution
simDiffusion=function(N,v,a,ter,z,s,dt,maxiter) {
  resp = rep(-1,N)
  rt = rep(-1,N)
  x = array(NA,c(N,maxiter+1))
  x[,1] = z
  for (n in 1:N) {
    for (i in 1:maxiter) {
      x[n,i+1] = x[n,i] + v*dt + s*sqrt(dt)*rnorm(1)
      if (x[n,i+1] > a) {
        resp[n] = 2
        rt[n] = ter + i*dt - dt/2
        break
      } 
      if (x[n,i+1] < 0) {
        resp[n] = 1
        rt[n] = ter + i*dt - dt/2
        break
      }
    }
  }
  return(list(resp=resp,rt=rt,x=x))
}


N=100
v=1
a=2
ter=0.3
z=a/2
s=1
dt=0.01
maxiter=1000

tmp=simDiffusion(N=N,v=v,a=a,ter=ter,z=z,s=s,dt=dt,maxiter=maxiter)

hist(tmp$rt)
plot(density(tmp$rt))
cat("Correct Proportion: ",mean(tmp$resp==2),"\n")
cat("Error Proportion: ",mean(tmp$resp==1),"\n")

plot(dt*1:(maxiter+1),tmp$x[1,],xlab="Time",ylab="Evidence",type="l",ylim=c(-a/2,a+a/2),col="grey",lty=2)
for (i in 2:N) {
  lines(dt*1:(maxiter+1),tmp$x[i,],col="grey",lty=2)
}
abline(h=a,col="red")
abline(h=0,col="red")


v=2
tmp=simDiffusion(N=N,v=v,a=a,ter=ter,z=z,s=s,dt=dt,maxiter=maxiter)
hist(tmp$rt)
plot(density(tmp$rt))
cat("Correct Proportion: ",mean(tmp$resp==2),"\n")
cat("Error Proportion: ",mean(tmp$resp==1),"\n")
plot(dt*1:(maxiter+1),tmp$x[1,],xlab="Time",ylab="Evidence",type="l",ylim=c(-a/2,a+a/2),col="grey",lty=2)
for (i in 2:N) {
  lines(dt*1:(maxiter+1),tmp$x[i,],col="grey",lty=2)
}
abline(h=a,col="red")
abline(h=0,col="red")


v=0.5
a=1
z=a/2
tmp=simDiffusion(N=N,v=v,a=a,ter=ter,z=z,s=s,dt=dt,maxiter=maxiter)

hist(tmp$rt)
plot(density(tmp$rt))
cat("Correct Proportion: ",mean(tmp$resp==2),"\n")
cat("Error Proportion: ",mean(tmp$resp==1),"\n")
plot(dt*1:(maxiter+1),tmp$x[1,],xlab="Time",ylab="Evidence",type="l",ylim=c(-a/2,a+a/2),col="grey",lty=2)
for (i in 2:N) {
  lines(dt*1:(maxiter+1),tmp$x[i,],col="grey",lty=2)
}
abline(h=a,col="red")
abline(h=0,col="red")


### step 5

rt1 <- rdiffusion(1000, a=1, v=2, t0=0.5, d=0, sz=0, sv=0, st0=0)
rt1 %>% group_by(response) %>% summarise(mrt = mean(rt))

rt2 <- rdiffusion(1000, a=2, v=2, t0=0.5, d=0, sz=0, sv=0, st0=0)
rt2 %>% group_by(response) %>% summarise(mrt = mean(rt))

rt3 <- rdiffusion(1000, a=2, v=3, t0=0.5, d=0, sz=0, sv=0, st0=0)
rt3 %>% group_by(response) %>% summarise(mrt = mean(rt))


# solution 

N=10000
ter=0.3


allV=seq(0,3,0.2)
allA=seq(0.2,3,0.2)

MCRT=array(NA,c(length(allV),length(allA)))
MERT=array(NA,c(length(allV),length(allA)))
ErrorProp=array(NA,c(length(allV),length(allA)))

for (vIndex in 1:length(allV)) {
  for (aIndex in 1:length(allA)) {
    v=allV[vIndex]
    a=allA[aIndex]
    tmp=rdiffusion(n=N,a=a,v=v,t0=ter)
    MCRT[vIndex,aIndex]=mean(tmp$rt[tmp$response=="upper"])
    MERT[vIndex,aIndex]=mean(tmp$rt[tmp$response=="lower"])
    ErrorProp[vIndex,aIndex]=mean(tmp$response=="lower")
  }
}

plot(MCRT,MERT,xlim=c(0.3,1.5),ylim=c(0.3,1.5),pch=16)
lines(x=c(0,100),y=c(0,100),col="red")

plot(allV,MCRT[,1]-MERT[,1],xlim=c(-0.5,3.5),ylim=c(-0.3,0.3),xlab="v",ylab="MRT Difference",pch=16)
for (i in 2:length(allA)) points(allV,MCRT[,i]-MERT[,i],pch=16)
abline(h=0,col="red")

plot(allA,MCRT[1,]-MERT[1,],xlim=c(0,3.2),ylim=c(-0.3,0.3),xlab="a",ylab="MRT Difference",pch=16)
for (i in 2:length(allV)) points(allA,MCRT[i,]-MERT[i,],pch=16)
abline(h=0,col="red")

plot(ErrorProp,MCRT-MERT,xlim=c(0,1),ylim=c(-0.3,0.3),xlab="Error Proportion",ylab="MRT Difference",pch=16)
abline(h=0,col="red")

# Session 2

logLikelihoodFunction=function(parameters, # parameters that should be fitted/optimized
                               parameterNames,
                               rt,
                               resp) {
  
  names(parameters) = parameterNames
  
  a = parameters["a"]
  v = parameters["v"]
  t0 = parameters["t0"]
  
  #Q: why logs
  #Q: why pmax and 1e-10
  out = sum(log(pmax( #pmax guards against -Inf likelihoods (if lh below 1e-10, set to 1e-10)
    ddiffusion(rt = rt, response = resp, # ddiffusion gives the likelihood for a certain parameter combination
               a = a, v = v, t0 = t0, 
               s = 0.1), 1e-10))) # s is a scaling parameter (diffusion constant)
  
  #Q: why negative?
  -out # optimizers usually minimize
  
}

#read in the data
dat = read.table(file = 'data/yo_dat.txt', header = T)

#here we add 1 because the objective function needs the responses to be either 'upper' and 'lower' or 1 and 2
responses = dat$choice + 1
dat = cbind(dat, responses)
hist(dat$rt)

#see how the likelihood function works
logLikelihoodFunction(parameters = c(0.3,0.2,0.2), 
                      parameterNames = c("a","v","t0"), 
                      rt = dat$rt,
                      resp = responses)

#Q: What are we going to do with this?

#starting parameters
startPar = c(0.2,0.2,0.2)
parameterNames = c('a','v','t0')
#optimise to find the best fitting parameters
# optim = simplex
outs = optim(par = startPar, logLikelihoodFunction, 
             parameterNames = parameterNames, 
             rt = dat$rt,
             resp = dat$responses)

fittedPars = c(outs$par[1], outs$par[2], outs$par[3])
options(scipen =999)
round(fittedPars, 4)

# fit data separately for conditions
outs = optim(par = startPar, logLikelihoodFunction, 
             parameterNames = parameterNames, 
             rt = dat$rt[dat$conds == 2], 
             resp = dat$responses)


#Q: What's dumb about what we just did?

#Break up the data by condition

#Break up the data by subject and fit all the individuals separately

#compare the older and younger individuals



# fit data separately for condition


N_sub <- length(unique(dat$subjs))
params_sub <- vector("list", N_sub)
params_sub
N_cond <- 2
params_cond <- vector("list", N_cond)

for(cond in 1:2){ 
  
  for(sub in seq_along(unique(dat$subjs))){ 
  
outs = optim(par = startPar, logLikelihoodFunction, 
             parameterNames = parameterNames, 
             rt = dat$rt[dat$conds == cond & dat$subjs == sub], 
             resp = dat$responses[dat$conds == cond & dat$subjs == sub])

params_sub[[sub]] = tibble(parameter = c("a", "v", "t0"), 
                       value = outs$par,
                       subject = sub, 
                       condition = cond)
  }
  params_sub_df <- bind_rows(params_sub) %>% pivot_wider(names_from = "parameter")
  params_cond[[cond]] <- params_sub_df
}
params_cond_df <- bind_rows(params_cond)
params_cond_df %>% View()

params_cond_df




## MCMC sampler  

# own solution 

# simulate data
dat <- rnorm(1000, 2, 2)
dat

# sampler 

mcmc <- function(data, iter, burnin, start, chains){ 
  
  posterior <- vector("numeric", iter) # vector for posteriors samples
  posterior[1] <- start # starting value
  
  for(i in 2:iter){ 
  
  prop_mean <- posterior[i-1] + rnorm(1, 0, .5)
  
  # check likelihood for proposal
  lh_prop <- sum(log(dnorm(data, mean = prop_mean, 2)))
  lh_last <- sum(log(dnorm(data, mean = posterior[i-1], 2)))
  
  # accept or not
  
  if(lh_prop > lh_last){ 
    posterior[i] <- prop_mean
    } else{
      prob_accept <-  exp(lh_prop - lh_last) 
      posterior[i] <- sample(c(prop_mean, posterior[i-1]), 1, prob = c(prob_accept, (1-prob_accept)))
      }
  }
return(posterior)
}

samples <- mcmc(data=dat, iter = 10000, burnin = 0, start = 2)
samples
hist(samples, breaks = 30)


# instructor solution 1 parameter

rm(list=ls())


# Step 1: A model with a likelihood function

logLikelihoodFunction=function(x,mean) {
  sum(log(dnorm(x,mean,fixedSD)))
}


# Step 2: Data

fixedSD=2
generatingMean=1

data=rnorm(1000,generatingMean,fixedSD)


# Step 3: Tuning parameters

iterations=5000
burnin=2000
chains=12

jumpSD=0.5

makeJump=function(n) {
  rnorm(n,0,jumpSD)
}


# Step 4: Storage

theta=array(NA,c(chains,iterations))
weight=rep(-Inf,chains)


# Step 5: Reasonable starting values

for (j in 1:chains) {
  while(weight[j]==-Inf) {
    theta[j,1]=rnorm(1,0,3)
    weight[j]=logLikelihoodFunction(data,theta[j,1])
  }
}


# Step 6: Iterative process

for (i in 2:iterations) {
  for (j in 1:chains) {
    oldParam=theta[j,i-1]
    oldWeight=weight[j]
    newParam=oldParam+makeJump(1)
    newWeight=logLikelihoodFunction(data,newParam)
    if(exp(newWeight-oldWeight) > runif(1)) {
      theta[j,i]=newParam
      weight[j]=newWeight
    } else {
      theta[j,i]=oldParam
      weight[j]=oldWeight
    }
  }
}


hist(theta[,burnin:iterations], breaks = 30)



# instructor solution 2 parameters 


rm(list=ls())



# Step 1: A model with a likelihood function

logLikelihoodFunction=function(x,mean,sd) {
  sum(log(dnorm(x,mean,sd)))
}

parameters=c("mean","sd")
npars=length(parameters)


# Step 2: Data

generatingMean=1
generatingSD=2

data=rnorm(1000,generatingMean,generatingSD)


# Step 3: Tuning parameters

iterations=5000
burnin=2000
chains=12

jumpSD=0.5

makeJump=function(n) {
  rnorm(n,0,jumpSD)
}


# Step 4: Storage

theta=array(NA,c(chains,npars,iterations))
colnames(theta)=parameters
weight=rep(-Inf,chains)


# Step 5: Reasonable starting values

for (j in 1:chains) {
  while(weight[j]==-Inf) {
    theta[j,"mean",1]=rnorm(1,0,3)
    theta[j,"sd",1]=rtnorm(1,0,3,0,Inf)
    weight[j]=logLikelihoodFunction(data,theta[j,"mean",1],theta[j,"sd",1])
  }
}


# Step 6: Iterative process

for (i in 2:iterations) {
  for (j in 1:chains) {
    oldParams=theta[j,,i-1]
    oldWeight=weight[j]
    newParams=oldParams+makeJump(npars)
    if (newParams[2]<0) {
      theta[j,,i]=oldParams
      weight[j]=oldWeight
    } else {
      newWeight=logLikelihoodFunction(data,newParams[1],newParams[2])
      if(exp(newWeight-oldWeight) > runif(1)) {
        theta[j,,i]=newParams
        weight[j]=newWeight
      } else {
        theta[j,,i]=oldParams
        weight[j]=oldWeight
      }
    }
  }
}


hist(theta[,1,burnin:iterations])
hist(theta[,2,burnin:iterations])


#step 3 priors

rm(list=ls())
library(msm)


# Step 1: A model with a likelihood function

logLikelihoodFunction=function(x,mean,sd) {
  sum(log(dnorm(x,mean,sd)))
}

logPriorMean=function(x){
  sum(log(dnorm(x,0,3)))
}

logPriorSD=function(x){
  sum(log(dgamma(x,1,1)))
}

parameters=c("mean","sd")
npars=length(parameters)


# Step 2: Data

generatingMean=1
generatingSD=2

data=rnorm(1000,generatingMean,generatingSD)


# Step 3: Tuning parameters

iterations=5000
burnin=2000
chains=12

jumpSD=0.5

makeJump=function(n) {
  rnorm(n,0,jumpSD)
}


# Step 4: Storage

theta=array(NA,c(chains,npars,iterations))
colnames(theta)=parameters
weight=rep(-Inf,chains)


# Step 5: Reasonable starting values

for (j in 1:chains) {
  while(weight[j]==-Inf) {
    theta[j,"mean",1]=rnorm(1,0,3)
    theta[j,"sd",1]=rtnorm(1,0,3,0,Inf)
    weight[j]=logLikelihoodFunction(data,theta[j,"mean",1],theta[j,"sd",1]) +
      logPriorMean(theta[j,"mean",1]) +
      logPriorSD(theta[j,"sd",1])
  }
}


# Step 6: Iterative process

for (i in 2:iterations) {
  for (j in 1:chains) {
    oldParams=theta[j,,i-1]
    oldWeight=weight[j]
    newParams=oldParams+makeJump(npars)
    if (newParams[2]<0) {
      theta[j,,i]=oldParams
      weight[j]=oldWeight
    } else {
      newWeight=logLikelihoodFunction(data,newParams[1],newParams[2]) +
        logPriorMean(newParams[1]) +
        logPriorSD(newParams[2])
      if(exp(newWeight-oldWeight) > runif(1)) {
        theta[j,,i]=newParams
        weight[j]=newWeight
      } else {
        theta[j,,i]=oldParams
        weight[j]=oldWeight
      }
    }
  }
}


hist(theta[,1,burnin:iterations])
hist(theta[,2,burnin:iterations])


# use differential evolution algorithm 

# own solution 


rm(list=ls())
library(msm)


# Step 1: A model with a likelihood function

logLikelihoodFunction=function(x,mean,sd) {
  sum(log(dnorm(x,mean,sd)))
}

logPriorMean=function(x){
  sum(log(dnorm(x,0,3)))
}

logPriorSD=function(x){
  sum(log(dgamma(x,1,1)))
}

parameters=c("mean","sd")
npars=length(parameters)


# Step 2: Data

generatingMean=1
generatingSD=2

data=rnorm(1000,generatingMean,generatingSD)


# Step 3: Tuning parameters

iterations=5000
burnin=2000
chains=12

# DE jumps
# jumpSD=0.5

# makeJump=function(n) {
  
#  rnorm(n,0,jumpSD)
  
#}


# Step 4: Storage

theta=array(NA,c(chains,npars,iterations))
colnames(theta)=parameters
weight=rep(-Inf,chains)


# Step 5: Reasonable starting values

for (j in 1:chains) {
  while(weight[j]==-Inf) {
    theta[j,"mean",1]=rnorm(1,0,3)
    theta[j,"sd",1]=rtnorm(1,0,3,0,Inf)
    weight[j]=logLikelihoodFunction(data,theta[j,"mean",1],theta[j,"sd",1]) +
      logPriorMean(theta[j,"mean",1]) +
      logPriorSD(theta[j,"sd",1])
  }
}


# Step 6: Iterative process

for (i in 2:iterations) {
  for (j in 1:chains) {
    oldParams=theta[j,,i-1]
    oldWeight=weight[j]
    newMean = oldParams[1] + diff(sample(theta[-j,"mean",i-1], 2)) # CAVE: Samples of mean and sd can be drawn from different pairs
    newSd = oldParams[2] + diff(sample(theta[-j,"sd",i-1], 2)) # CAVE: Samples of mean and sd can be drawn from different pairs
    newParams = c(newMean, newSd)
    if (newParams[2]<0) {
      theta[j,,i]=oldParams
      weight[j]=oldWeight
    } else {
      newWeight=logLikelihoodFunction(data,newParams[1],newParams[2]) +
        logPriorMean(newParams[1]) +
        logPriorSD(newParams[2])
      if(exp(newWeight-oldWeight) > runif(1)) {
        theta[j,,i]=newParams
        weight[j]=newWeight
      } else {
        theta[j,,i]=oldParams
        weight[j]=oldWeight
      }
    }
  }
}


hist(theta[,1,burnin:iterations])
hist(theta[,2,burnin:iterations])

# instructor solution 

m(list=ls())
library(msm)


# Step 1: A model with a likelihood function

logLikelihoodFunction=function(x,mean,sd) {
  sum(log(dnorm(x,mean,sd)))
}

logPriorMean=function(x){
  sum(log(dnorm(x,0,3)))
}

logPriorSD=function(x){
  sum(log(dgamma(x,1,1)))
}

parameters=c("mean","sd")
npars=length(parameters)


# Step 2: Data

generatingMean=1
generatingSD=2

data=rnorm(1000,generatingMean,generatingSD)


# Step 3: Tuning parameters

iterations=5000
burnin=2000
chains=12

jumpSD=0.5

makeJump=function(otherChainValues1,otherChainValues2) {
  otherChainValues1-otherChainValues2
}


# Step 4: Storage

theta=array(NA,c(chains,npars,iterations))
colnames(theta)=parameters
weight=rep(-Inf,chains)


# Step 5: Reasonable starting values

for (j in 1:chains) {
  while(weight[j]==-Inf) {
    theta[j,"mean",1]=rnorm(1,0,3)
    theta[j,"sd",1]=rtnorm(1,0,3,0,Inf)
    weight[j]=logLikelihoodFunction(data,theta[j,"mean",1],theta[j,"sd",1]) +
      logPriorMean(theta[j,"mean",1]) +
      logPriorSD(theta[j,"sd",1])
  }
}


# Step 6: Iterative process

for (i in 2:iterations) {
  for (j in 1:chains) {
    oldParams=theta[j,,i-1]
    oldWeight=weight[j]
    allChains=1:chains
    otherChains=sample(allChains[-j],2)
    otherChain1=otherChains[1]
    otherChain2=otherChains[2]
    newParams=oldParams+makeJump(theta[otherChain1,,i-1],theta[otherChain2,,i-1])
    if (newParams[2]<0) {
      theta[j,,i]=oldParams
      weight[j]=oldWeight
    } else {
      newWeight=logLikelihoodFunction(data,newParams[1],newParams[2]) +
        logPriorMean(newParams[1]) +
        logPriorSD(newParams[2])
      if(exp(newWeight-oldWeight) > runif(1)) {
        theta[j,,i]=newParams
        weight[j]=newWeight
      } else {
        theta[j,,i]=oldParams
        weight[j]=oldWeight
      }
    }
  }
}


hist(theta[,1,burnin:iterations])
hist(theta[,2,burnin:iterations])


# MCMC for diffusion model 


# own solution 

m(list=ls())
library(msm)


# Step 1: A model with a likelihood function

logLikelihoodFunction=function(parameters, parameterNames, rt, resp) {
  
  names(parameters) = parameterNames
  
  a = parameters["a"]
  v = parameters["v"]
  t0 = parameters["t0"]
  
  out = sum(log(pmax(
    ddiffusion(rt = rt, response = resp, a = a, v = v, t0 = t0, s = 0.1), 1e-10)))
  
}

logPrior_v=function(x){
  sum(log(dnorm(x,0,2.5)))
}

logPrior_a=function(x){
  sum(log(dgamma(x,1,1)))
}

logPrior_t0=function(x){ 
  sum(log(dgamma(x, 2, 2)))
  }

parameters=c("a","v", "t0")
npars=length(parameters)


# Step 2: Data

generating_a=2
generating_v=1
generating_t0=.3

data=rdiffusion(1000,t0 = generating_t0, a = generating_a, v = generating_v)


# Step 3: Tuning parameters

iterations=5000
burnin=2000
chains=12

jumpSD=0.5

makeJump=function(otherChainValues1,otherChainValues2) {
  otherChainValues1-otherChainValues2
}


# Step 4: Storage

theta=array(NA,c(chains,npars,iterations))
colnames(theta)=parameters
weight=rep(-Inf,chains)


# Step 5: Reasonable starting values

for (j in 1:chains) {
  while(weight[j]==-Inf) {
    theta[j,"mean",1]=rnorm(1,0,3)
    theta[j,"sd",1]=rtnorm(1,0,3,0,Inf)
    weight[j]=logLikelihoodFunction(rt = data$rt, 
                                    response = data$response, 
                                    theta[j,"mean",1],theta[j,"sd",1]) +
      logPriorMean(theta[j,"mean",1]) +
      logPriorSD(theta[j,"sd",1])
  }
}


# Step 6: Iterative process

for (i in 2:iterations) {
  for (j in 1:chains) {
    oldParams=theta[j,,i-1]
    oldWeight=weight[j]
    allChains=1:chains
    otherChains=sample(allChains[-j],2)
    otherChain1=otherChains[1]
    otherChain2=otherChains[2]
    newParams=oldParams+makeJump(theta[otherChain1,,i-1],theta[otherChain2,,i-1])
    if (newParams[2]<0) {
      theta[j,,i]=oldParams
      weight[j]=oldWeight
    } else {
      newWeight=logLikelihoodFunction(data,newParams[1],newParams[2]) +
        logPriorMean(newParams[1]) +
        logPriorSD(newParams[2])
      if(exp(newWeight-oldWeight) > runif(1)) {
        theta[j,,i]=newParams
        weight[j]=newWeight
      } else {
        theta[j,,i]=oldParams
        weight[j]=oldWeight
      }
    }
  }
}


hist(theta[,1,burnin:iterations])
hist(theta[,2,burnin:iterations])




# instructor solution 

rm(list=ls())
library(msm)

# Step 1: A model with a likelihood function

logLikelihoodFunction=function(parameters,parameterNames,rt,resp) {
  
  names(parameters) = parameterNames
  out=0
  a = parameters["a"]
  v = parameters["v"]
  t0 = parameters["t0"]
  out = out + sum(log(pmax(
    ddiffusion(rt = rt, response = resp,
               a = a, v = v, t0 = t0),1e-10)))
  out
}

logPriorV=function(x){
  sum(log(dtnorm(x,3,3,0,Inf)))
}

logPriorA=function(x){
  sum(log(dtnorm(x,2,2,0,Inf)))
}

logPriorT0=function(x){
  sum(log(dtnorm(x,0.3,0.3,0,Inf)))
}

parameterNames=c("v","a","t0")
npars=length(parameterNames)


# Step 2: Data

generatingV=1
generatingA=2
generatingT0=0.3


responseTimes=NULL
responses=NULL

useN=1000


data=rdiffusion(n=useN, a = generatingA, v = generatingV, t0 = generatingT0)

responseTimes=c(responseTimes,as.numeric(data$rt))
responses=c(responses,as.numeric(data$response))


# Step 3: Tuning parameters

iterations=2000
burnin=500
chains=12

makeJump=function(prevParams,otherChainValues1,otherChainValues2) {
  out = NA
  useNPars = length(prevParams)
  gamma = 2.38/sqrt(2*useNPars)
  mutation = 0.0001
  out = prevParams + gamma*(otherChainValues1-otherChainValues2) + runif(useNPars,-mutation,mutation)
  out
}


# Step 4: Storage

theta=array(NA,c(chains,npars,iterations))
colnames(theta)=parameterNames
weight=rep(-Inf,chains)
storeMLE=array(NA,c(chains,iterations))


# Step 5: Reasonable starting values

for (j in 1:chains) {
  while(weight[j]==-Inf) {
    theta[j,"v",1] = rtnorm(1,3,3,0,Inf)
    theta[j,"a",1] = rtnorm(1,2,2,0,Inf)
    theta[j,"t0",1] = rtnorm(1,0.3,0.3,0,Inf)
    tmpLike = logLikelihoodFunction(theta[j,,1],parameterNames,
                                    responseTimes,responses)
    
    weight[j] = tmpLike + 
      logPriorV(theta[j,"v",1]) +
      logPriorA(theta[j,"a",1]) +
      logPriorT0(theta[j,"t0",1])
  }
  storeMLE[j,1] = tmpLike
}


# Step 6: Iterative process

for (i in 2:iterations) {
  cat(i," ")
  for (j in 1:chains) {
    oldParams=theta[j,,i-1]
    oldWeight=weight[j]
    allChains=1:chains
    otherChains=sample(allChains[-j],2)
    otherChain1=otherChains[1]
    otherChain2=otherChains[2]
    newParams=makeJump(oldParams,theta[otherChain1,,i-1],theta[otherChain2,,i-1])
    names(newParams)=parameterNames
    if (any(newParams<0)) {
      theta[j,,i]=oldParams
      weight[j]=oldWeight
      storeMLE[j,i]=storeMLE[j,i-1]
    } else {
      tmpLike = logLikelihoodFunction(newParams,parameterNames,
                                      responseTimes,responses)
      
      newWeight=tmpLike + 
        logPriorV(newParams["v"]) +
        logPriorA(newParams["a"]) +
        logPriorT0(newParams["t0"])
      if(exp(newWeight-oldWeight) > runif(1)) {
        theta[j,,i]=newParams
        weight[j]=newWeight
        storeMLE[j,i]=tmpLike
      } else {
        theta[j,,i]=oldParams
        weight[j]=oldWeight
        storeMLE[j,i]=storeMLE[j,i-1]
      }
    }
  }
}

hist(theta[,"v",burnin:iterations])
hist(theta[,"a",burnin:iterations])
hist(theta[,"t0",burnin:iterations])


