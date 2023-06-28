rm(list = ls())
logLikelihoodFunction=function(parameters,parameterNames,
                               rt,resp,conditions) {
  
  #what does this new variable 'conditions' do?
  nconditions = length(unique(conditions))
  
  names(parameters) = parameterNames
  out=0
  
  #what is going on here?
  a = parameters[grep("a",parameterNames)]
  if (length(a) < nconditions) a = rep(a, nconditions)
  v = parameters[grep("v",parameterNames)]
  if (length(v) < nconditions) v = rep(v, nconditions)
  t0 = parameters[grep("t0",parameterNames)]
  if (length(t0) < nconditions) t0 = rep(t0, nconditions)
  
  #notice the loop here
  for (i in 1:nconditions){
    out = out + sum(log(pmax(
      ddiffusion(rt = rt[conditions == i], response = resp[conditions == i],
                 a = a[i], v = v[i], t0 = t0[i], 
                 s = 0.1), 1e-10)))
  }
  -out
  
}

dat = read.table(file = 'data/yo_dat_conds.txt', header = T)
responses = dat$choice + 1
dat = cbind(dat, responses)

ngroups = length(unique(dat$group))
nsubjs = length(unique(dat$subjs))

#more parameters per participant - you should understand why
startPar = c(0.1, 0.1, 0.2, 0.3, 0.15)
#why are these parameter names important?
parameterNames = c("a", "v1", "v2", "v3", "t0")
####how would you change this to let a different parameter vary across conditions?
fittedPars = array(NA, dim = c(5, nsubjs, ngroups))
for (i in 1:ngroups){
  for (j in 1:nsubjs){
    print(j)
    #where does the conditions variable play its part?
    tmpdat = dat[dat$group == i & dat$subjs == j & !is.na(dat$rt), ]
    outs = optim(par = startPar, logLikelihoodFunction, 
                 parameterNames = parameterNames, 
                 rt = tmpdat$rt, 
                 resp = tmpdat$responses,
                 conditions = tmpdat$conds)
    
    fittedPars[,j,i] = outs$par
  }
}

mns = apply(fittedPars, c(1,3), mean)
sds = apply(fittedPars, c(1,3), sd)

layout(m = array(1:5, dim = c(1,5)))
par(mar = c(4,4,3,1))
plot(mns[1,], ylim = c(0.1,0.4), xlim = c(0.75, 2.25), axes = F, xlab = '', ylab = '', main = 'boundary separation')
axis(side = 1, at = 1:2, c('young','old')); axis(side = 2); box()
arrows(x0 = 1:2, y0 = mns[1,]-sds[1,]/sqrt(nsubjs), y1 = mns[1,] + sds[1,]/sqrt(nsubjs), code = 3, angle = 90, length = 0.05)
plot(mns[2,], ylim = c(0,0.4), xlim = c(0.75, 2.25), axes = F, xlab = '', ylab = '', main = 'drift rate')
axis(side = 1, at = 1:2, c('young','old')); axis(side = 2); box()
arrows(x0 = 1:2, y0 = mns[2,]-sds[2,]/sqrt(nsubjs), y1 = mns[2,] + sds[2,]/sqrt(nsubjs), code = 3, angle = 90, length = 0.05)
plot(mns[3,], ylim = c(0,0.4), xlim = c(0.75, 2.25), axes = F, xlab = '', ylab = '', main = 'drift rate')
axis(side = 1, at = 1:2, c('young','old')); axis(side = 2); box()
arrows(x0 = 1:2, y0 = mns[3,]-sds[3,]/sqrt(nsubjs), y1 = mns[3,] + sds[3,]/sqrt(nsubjs), code = 3, angle = 90, length = 0.05)
plot(mns[4,], ylim = c(0,0.4), xlim = c(0.75, 2.25), axes = F, xlab = '', ylab = '', main = 'drift rate')
axis(side = 1, at = 1:2, c('young','old')); axis(side = 2); box()
arrows(x0 = 1:2, y0 = mns[4,]-sds[4,]/sqrt(nsubjs), y1 = mns[4,] + sds[4,]/sqrt(nsubjs), code = 3, angle = 90, length = 0.05)
plot(mns[5,], ylim = c(0,0.35), xlim = c(0.75, 2.25), axes = F, xlab = '', ylab = '', main = 'non-decision time')
axis(side = 1, at = 1:2, c('young','old')); axis(side = 2); box()
arrows(x0 = 1:2, y0 = mns[5,]-sds[5,]/sqrt(nsubjs), y1 = mns[5,] + sds[5,]/sqrt(nsubjs), code = 3, angle = 90, length = 0.05)

t.test(fittedPars[1,,1], fittedPars[1,,2])
t.test(fittedPars[2,,1], fittedPars[2,,2])
t.test(fittedPars[3,,1], fittedPars[3,,2])
t.test(fittedPars[4,,1], fittedPars[4,,2])
t.test(fittedPars[5,,1], fittedPars[5,,2])



# Mystery Data Set 1 ------------------------------------------------------

# Exercise 1a: Fit model simultaenously for both conditions 


rm(list=ls())
library(msm)

useDataset=1

# Step 1: A model with a likelihood function

logLikelihoodFunction=function(parameters,parameterNames,rt,resp,
                               cond=NULL,conditions=NULL) {
  
  names(parameters) = parameterNames
  out=0
  if (is.null(conditions)) {
    a = parameters["a"]
    v = parameters["v"]
    t0 = parameters["t0"]
    out = out + sum(log(pmax(
      ddiffusion(rt = rt, response = resp,
                 a = a, v = v, t0 = t0),1e-10)))
  } else {
    for (useCond in conditions) {
      isCond = cond==useCond
      a = parameters["a"]
      v = parameters["v"]
      t0 = parameters["t0"]
      out = out + sum(log(pmax(
        ddiffusion(rt = rt[isCond], response = resp[isCond],
                   a = a, v = v, t0 = t0),1e-10)))
    }
  }
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

load(paste("data/exercise1_mysteryDataset",useDataset,".Rdata",sep=""))


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
                                    responseTimes,responses,
                                    cond=cond,conditions=conditions)
    
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
                                      responseTimes,responses,
                                      cond=cond,conditions=conditions)
      
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

AIC = max(storeMLE)*(-2) + 2*npars
BIC = max(storeMLE)*(-2) + log(length(responseTimes))*npars
aveDev = mean(storeMLE[,burnin:iterations])*(-2)
maxDev = max(storeMLE[,burnin:iterations])*(-2)
DIC = aveDev + (aveDev-maxDev)

save(AIC,BIC,DIC,file=paste("modelSelection_exercise1_withData",useDataset,"_sim.Rdata",sep=""))
M1 <- data.frame(AIC, BIC, DIC, M = "1")
M1
# Exercise 1b: Thresholds are allowed to vary

rm(list=ls())
library(msm)

useDataset=1

# Step 1: A model with a likelihood function

logLikelihoodFunction=function(parameters,parameterNames,rt,resp,
                               cond=NULL,conditions=NULL) {
  
  names(parameters) = parameterNames
  out=0
  if (is.null(conditions)) {
    a = parameters["a"]
    v = parameters["v"]
    t0 = parameters["t0"]
    out = out + sum(log(pmax(
      ddiffusion(rt = rt, response = resp,
                 a = a, v = v, t0 = t0),1e-10)))
  } else {
    for (useCond in conditions) {
      isCond = cond==useCond
      a = parameters[paste("a",useCond,sep=".")]
      v = parameters["v"]
      t0 = parameters["t0"]
      out = out + sum(log(pmax(
        ddiffusion(rt = rt[isCond], response = resp[isCond],
                   a = a, v = v, t0 = t0),1e-10)))
    }
  }
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

parameterNames=c("v","a.accuracy","a.speed","t0")
npars=length(parameterNames)

# Step 2: Data

load(paste("data/exercise1_mysteryDataset",useDataset,".Rdata",sep=""))


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
    theta[j,"a.accuracy",1] = rtnorm(1,2,2,0,Inf)
    theta[j,"a.speed",1] = rtnorm(1,2,2,0,Inf)
    theta[j,"t0",1] = rtnorm(1,0.3,0.3,0,Inf)
    tmpLike = logLikelihoodFunction(theta[j,,1],parameterNames,
                                    responseTimes,responses,
                                    cond=cond,conditions=conditions)
    
    weight[j] = tmpLike + 
      logPriorV(theta[j,"v",1]) +
      logPriorA(theta[j,"a.accuracy",1]) +
      logPriorA(theta[j,"a.speed",1]) +
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
                                      responseTimes,responses,
                                      cond=cond,conditions=conditions)
      
      newWeight=tmpLike + 
        logPriorV(newParams["v"]) +
        logPriorA(newParams["a.accuracy"]) +
        logPriorA(newParams["a.speed"]) +
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
hist(theta[,"a.accuracy",burnin:iterations])
hist(theta[,"a.speed",burnin:iterations])
hist(theta[,"t0",burnin:iterations])

AIC = max(storeMLE)*(-2) + 2*npars
BIC = max(storeMLE)*(-2) + log(length(responseTimes))*npars
aveDev = mean(storeMLE[,burnin:iterations])*(-2)
maxDev = max(storeMLE[,burnin:iterations])*(-2)
DIC = aveDev + (aveDev-maxDev)

save(AIC,BIC,DIC,file=paste("modelSelection_exercise1_withData",useDataset,"_avary.Rdata",sep=""))
M2 <- data.frame(AIC, BIC, DIC, M = "2")
M2
# Exercise 1c: All parameters are allowed to vary

library(msm)

useDataset=1

# Step 1: A model with a likelihood function

logLikelihoodFunction=function(parameters,parameterNames,rt,resp,
                               cond=NULL,conditions=NULL) {
  
  names(parameters) = parameterNames
  out=0
  if (is.null(conditions)) {
    a = parameters["a"]
    v = parameters["v"]
    t0 = parameters["t0"]
    out = out + sum(log(pmax(
      ddiffusion(rt = rt, response = resp,
                 a = a, v = v, t0 = t0),1e-10)))
  } else {
    for (useCond in conditions) {
      isCond = cond==useCond
      a = parameters[paste("a",useCond,sep=".")]
      v = parameters[paste("v",useCond,sep=".")]
      t0 = parameters["t0"]
      out = out + sum(log(pmax(
        ddiffusion(rt = rt[isCond], response = resp[isCond],
                   a = a, v = v, t0 = t0),1e-10)))
    }
  }
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

parameterNames=c("v.accuracy","v.speed","a.accuracy","a.speed","t0")
npars=length(parameterNames)

# Step 2: Data

load(paste("data/exercise1_mysteryDataset",useDataset,".Rdata",sep=""))


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
    theta[j,"v.accuracy",1] = rtnorm(1,3,3,0,Inf)
    theta[j,"v.speed",1] = rtnorm(1,3,3,0,Inf)
    theta[j,"a.accuracy",1] = rtnorm(1,2,2,0,Inf)
    theta[j,"a.speed",1] = rtnorm(1,2,2,0,Inf)
    theta[j,"t0",1] = rtnorm(1,0.3,0.3,0,Inf)
    tmpLike = logLikelihoodFunction(theta[j,,1],parameterNames,
                                    responseTimes,responses,
                                    cond=cond,conditions=conditions)
    
    weight[j] = tmpLike + 
      logPriorV(theta[j,"v.accuracy",1]) +
      logPriorV(theta[j,"v.speed",1]) +
      logPriorA(theta[j,"a.accuracy",1]) +
      logPriorA(theta[j,"a.speed",1]) +
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
                                      responseTimes,responses,
                                      cond=cond,conditions=conditions)
      
      newWeight=tmpLike + 
        logPriorV(newParams["v.accuracy"]) +
        logPriorV(newParams["v.speed"]) +
        logPriorA(newParams["a.accuracy"]) +
        logPriorA(newParams["a.speed"]) +
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


hist(theta[,"v.accuracy",burnin:iterations])
hist(theta[,"v.speed",burnin:iterations])
hist(theta[,"a.accuracy",burnin:iterations])
hist(theta[,"a.speed",burnin:iterations])
hist(theta[,"t0",burnin:iterations])

AIC = max(storeMLE)*(-2) + 2*npars
BIC = max(storeMLE)*(-2) + log(length(responseTimes))*npars
aveDev = mean(storeMLE[,burnin:iterations])*(-2)
maxDev = max(storeMLE[,burnin:iterations])*(-2)
DIC = aveDev + (aveDev-maxDev)

M3_D1 <- data.frame(AIC, BIC, DIC, M = "3")

D1 <- bind_rows(M1_D1, M2_D1, M3_D1)



# Mystery Data Set 2 ------------------------------------------------------

# Exercise 1a: Fit model simultaenously for both conditions 




useDataset=2

# Step 1: A model with a likelihood function

logLikelihoodFunction=function(parameters,parameterNames,rt,resp,
                               cond=NULL,conditions=NULL) {
  
  names(parameters) = parameterNames
  out=0
  if (is.null(conditions)) {
    a = parameters["a"]
    v = parameters["v"]
    t0 = parameters["t0"]
    out = out + sum(log(pmax(
      ddiffusion(rt = rt, response = resp,
                 a = a, v = v, t0 = t0),1e-10)))
  } else {
    for (useCond in conditions) {
      isCond = cond==useCond
      a = parameters["a"]
      v = parameters["v"]
      t0 = parameters["t0"]
      out = out + sum(log(pmax(
        ddiffusion(rt = rt[isCond], response = resp[isCond],
                   a = a, v = v, t0 = t0),1e-10)))
    }
  }
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

load(paste("data/exercise1_mysteryDataset",useDataset,".Rdata",sep=""))


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
                                    responseTimes,responses,
                                    cond=cond,conditions=conditions)
    
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
                                      responseTimes,responses,
                                      cond=cond,conditions=conditions)
      
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

AIC = max(storeMLE)*(-2) + 2*npars
BIC = max(storeMLE)*(-2) + log(length(responseTimes))*npars
aveDev = mean(storeMLE[,burnin:iterations])*(-2)
maxDev = max(storeMLE[,burnin:iterations])*(-2)
DIC = aveDev + (aveDev-maxDev)

save(AIC,BIC,DIC,file=paste("modelSelection_exercise1_withData",useDataset,"_sim.Rdata",sep=""))
M1 <- data.frame(AIC, BIC, DIC, M = "1")
M1
# Exercise 1b: Thresholds are allowed to vary

useDataset=2

# Step 1: A model with a likelihood function

logLikelihoodFunction=function(parameters,parameterNames,rt,resp,
                               cond=NULL,conditions=NULL) {
  
  names(parameters) = parameterNames
  out=0
  if (is.null(conditions)) {
    a = parameters["a"]
    v = parameters["v"]
    t0 = parameters["t0"]
    out = out + sum(log(pmax(
      ddiffusion(rt = rt, response = resp,
                 a = a, v = v, t0 = t0),1e-10)))
  } else {
    for (useCond in conditions) {
      isCond = cond==useCond
      a = parameters[paste("a",useCond,sep=".")]
      v = parameters["v"]
      t0 = parameters["t0"]
      out = out + sum(log(pmax(
        ddiffusion(rt = rt[isCond], response = resp[isCond],
                   a = a, v = v, t0 = t0),1e-10)))
    }
  }
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

parameterNames=c("v","a.accuracy","a.speed","t0")
npars=length(parameterNames)

# Step 2: Data

load(paste("data/exercise1_mysteryDataset",useDataset,".Rdata",sep=""))


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
    theta[j,"a.accuracy",1] = rtnorm(1,2,2,0,Inf)
    theta[j,"a.speed",1] = rtnorm(1,2,2,0,Inf)
    theta[j,"t0",1] = rtnorm(1,0.3,0.3,0,Inf)
    tmpLike = logLikelihoodFunction(theta[j,,1],parameterNames,
                                    responseTimes,responses,
                                    cond=cond,conditions=conditions)
    
    weight[j] = tmpLike + 
      logPriorV(theta[j,"v",1]) +
      logPriorA(theta[j,"a.accuracy",1]) +
      logPriorA(theta[j,"a.speed",1]) +
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
                                      responseTimes,responses,
                                      cond=cond,conditions=conditions)
      
      newWeight=tmpLike + 
        logPriorV(newParams["v"]) +
        logPriorA(newParams["a.accuracy"]) +
        logPriorA(newParams["a.speed"]) +
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
hist(theta[,"a.accuracy",burnin:iterations])
hist(theta[,"a.speed",burnin:iterations])
hist(theta[,"t0",burnin:iterations])

AIC = max(storeMLE)*(-2) + 2*npars
BIC = max(storeMLE)*(-2) + log(length(responseTimes))*npars
aveDev = mean(storeMLE[,burnin:iterations])*(-2)
maxDev = max(storeMLE[,burnin:iterations])*(-2)
DIC = aveDev + (aveDev-maxDev)

save(AIC,BIC,DIC,file=paste("modelSelection_exercise1_withData",useDataset,"_avary.Rdata",sep=""))
M2 <- data.frame(AIC, BIC, DIC, M = "2")
M2
# Exercise 1c: All parameters are allowed to vary

library(msm)

useDataset=2

# Step 1: A model with a likelihood function

logLikelihoodFunction=function(parameters,parameterNames,rt,resp,
                               cond=NULL,conditions=NULL) {
  
  names(parameters) = parameterNames
  out=0
  if (is.null(conditions)) {
    a = parameters["a"]
    v = parameters["v"]
    t0 = parameters["t0"]
    out = out + sum(log(pmax(
      ddiffusion(rt = rt, response = resp,
                 a = a, v = v, t0 = t0),1e-10)))
  } else {
    for (useCond in conditions) {
      isCond = cond==useCond
      a = parameters[paste("a",useCond,sep=".")]
      v = parameters[paste("v",useCond,sep=".")]
      t0 = parameters["t0"]
      out = out + sum(log(pmax(
        ddiffusion(rt = rt[isCond], response = resp[isCond],
                   a = a, v = v, t0 = t0),1e-10)))
    }
  }
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

parameterNames=c("v.accuracy","v.speed","a.accuracy","a.speed","t0")
npars=length(parameterNames)

# Step 2: Data

load(paste("data/exercise1_mysteryDataset",useDataset,".Rdata",sep=""))


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
    theta[j,"v.accuracy",1] = rtnorm(1,3,3,0,Inf)
    theta[j,"v.speed",1] = rtnorm(1,3,3,0,Inf)
    theta[j,"a.accuracy",1] = rtnorm(1,2,2,0,Inf)
    theta[j,"a.speed",1] = rtnorm(1,2,2,0,Inf)
    theta[j,"t0",1] = rtnorm(1,0.3,0.3,0,Inf)
    tmpLike = logLikelihoodFunction(theta[j,,1],parameterNames,
                                    responseTimes,responses,
                                    cond=cond,conditions=conditions)
    
    weight[j] = tmpLike + 
      logPriorV(theta[j,"v.accuracy",1]) +
      logPriorV(theta[j,"v.speed",1]) +
      logPriorA(theta[j,"a.accuracy",1]) +
      logPriorA(theta[j,"a.speed",1]) +
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
                                      responseTimes,responses,
                                      cond=cond,conditions=conditions)
      
      newWeight=tmpLike + 
        logPriorV(newParams["v.accuracy"]) +
        logPriorV(newParams["v.speed"]) +
        logPriorA(newParams["a.accuracy"]) +
        logPriorA(newParams["a.speed"]) +
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


hist(theta[,"v.accuracy",burnin:iterations])
hist(theta[,"v.speed",burnin:iterations])
hist(theta[,"a.accuracy",burnin:iterations])
hist(theta[,"a.speed",burnin:iterations])
hist(theta[,"t0",burnin:iterations])

AIC = max(storeMLE)*(-2) + 2*npars
BIC = max(storeMLE)*(-2) + log(length(responseTimes))*npars
aveDev = mean(storeMLE[,burnin:iterations])*(-2)
maxDev = max(storeMLE[,burnin:iterations])*(-2)
DIC = aveDev + (aveDev-maxDev)

M3 <- data.frame(AIC, BIC, DIC, M = "3")

data2 <- bind_rows(M1, M2, M3)
data2











