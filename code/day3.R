# sample from prior 
N <- 1e3
smp_a <- runif(N, .05, .35)
smp_v <- runif(N, 0, .5)
smp_t <- runif(N, .2, .6)

parameters <- data.frame(smp_a, smp_v, smp_t)
parameters

parameterNames=c("v","a","t0")


# define likelihood


logLikelihoodFunction=function(parameters, 
                               parameterNames,
                               rt,
                               resp) {
  
  names(parameters) = parameterNames
  
  a = parameters["a"]
  v = parameters["v"]
  t0 = parameters["t0"]
  
  out = sum(log(pmax(
    ddiffusion(rt = rt, response = resp,
               a = a, v = v, t0 = t0, 
               s = 0.1), 1e-10)))
  
  -out
  
}

# read data

dat = read.table(file = 'data/yo_dat.txt', header = T)
responses = dat$choice + 1
dat = cbind(dat, responses)

# calculate the likelihood

lh <- vector("numeric", length = N)

for(comb in seq_len(nrow(parameters))){
  
  for(j in 1:2){
    
    for(i in 1:length(unique(dat$subjs))){
      
      lh[comb] <- logLikelihoodFunction(parameters= c(parameters[comb, "smp_a"], parameters[comb, "smp_v"], parameters[comb, "smp_t"]),
                      parameterNames = parameterNames, 
                      rt = dat$rt[dat$conds == j & dat$subjs == i], 
                      resp = dat$responses[dat$conds == j & dat$subjs == i])
      
    }
  }
}

mean(exp(-lh))

###

#####.  Let's update the old simulation code we had #####

simDiffusion=function(N, v_mean, v_sd,a,ter,z,s,dt,maxiter) {
  resp = rep(-1,N)
  rt = rep(-1,N)
  x = array(NA,c(N,maxiter+1))
  x[,1] = z
  iter = 0
  v = rnorm(N, v_mean, v_sd)
  keep = rep(TRUE,N)
  while (iter < maxiter) {
    iter = iter + 1
    x[keep,iter+1] = x[keep,iter] + v[keep]*dt + s*sqrt(dt)*rnorm(sum(keep))
    if (any(x[keep,iter+1] > a)) {
      resp[keep & x[,iter+1]>a] = 2
      rt[keep & x[,iter+1]>a] = ter + iter*dt - dt/2
      keep[keep & x[,iter+1]>a] = FALSE
    } 
    if (any(x[keep,iter+1] < 0)) {
      resp[keep & x[,iter+1]<0] = 1
      rt[keep & x[,iter+1]<0] = ter + iter*dt - dt/2
      keep[keep & x[,iter+1]<0] = FALSE
    }
    if (sum(keep)==0) {
      break
    }
  }
  return(list(resp=resp,rt=rt,x=x))
}


N=1000
v_mean = 0
v_sd = .3
a=2
ter=0.3
z=a/2
s=1
dt=0.001
maxiter=10000

tmp1=simDiffusion(N=N,v_mean=v_mean, v_sd=v_sd,a=a,ter=ter,z=z,s=s,dt=dt,maxiter=maxiter)


cat("Correct Proportion: ",mean(tmp1$resp==2),"\n")
cat("Error Proportion: ",mean(tmp1$resp==1),"\n")
cat("Mean Correct RT: ",mean(tmp1$rt[tmp1$resp==2]),"\n")
cat("Mean Error RT: ",mean(tmp1$rt[tmp1$resp==1]),"\n")

############################################
##### Use rdiffusion to do this faster #####

N=10000
ter=0.3
sv = 1.8
sz = .5

allV=seq(0,3,0.2)
allA=seq(0.2,3,0.2)

MCRT=array(NA,c(length(allV),length(allA)))
MERT=array(NA,c(length(allV),length(allA)))
ErrorProp=array(NA,c(length(allV),length(allA)))

# inter-trial variability for drift rate
for (vIndex in 1:length(allV)) {
  for (aIndex in 1:length(allA)) {
    v=allV[vIndex]
    a=allA[aIndex]
    tmp=rdiffusion(n=N,a=a,v=v,t0=ter, sv = sv )
    MCRT[vIndex,aIndex]=mean(tmp$rt[tmp$response=="upper"])
    MERT[vIndex,aIndex]=mean(tmp$rt[tmp$response=="lower"])
    ErrorProp[vIndex,aIndex]=mean(tmp$response=="lower")
  }
}


# inter-trial variability for drift rate and starting point

N=10000
ter=0.3
sv = 1.8
sz = .5

V = 1.8
A = 1

for (vIndex in 1:length(allV)) {
  for (aIndex in 1:length(allA)) {
    v=V
    a=A
    tmp=rdiffusion(n=N,a=a,v=v,t0=ter, sv = sv, sz = sz )
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


# qualitative predictions 

logLikelihoodFunction=function(parameters,parameterNames,
                               rt,resp,conditions) {
  
  nconditions = length(unique(conditions))
  
  names(parameters) = parameterNames
  out=0
  
  a = parameters[grep("a",parameterNames)]
  if (length(a) < nconditions) a = rep(a, nconditions)
  v = parameters[grep("v",parameterNames)]
  if (length(v) < nconditions) v = rep(v, nconditions)
  t0 = parameters[grep("t0",parameterNames)]
  if (length(t0) < nconditions) t0 = rep(t0, nconditions)
  
  for (i in 1:nconditions){
    out = out + sum(log(pmax(
      ddiffusion(rt = rt[conditions == i], response = resp[conditions == i],
                 a = a, v = v[i], t0 = t0, 
                 s = 0.1), 1e-10)))
  }
  -out
  
}

#Read in the new data set
dat = read.table(file = 'data/dat_conds.txt', header = T)
responses = dat$choice + 1
dat = cbind(dat, responses)
nsubjs = length(unique(dat$subjs))

#Using code from other sessions, fit this data set letting drift rate vary across conditions
startPar = c(0.1, 0.1, 0.2, 0.3, 0.15)
parameterNames = c("a", "v1", "v2", "v3", "t0")

####how would you change this to let a different parameter vary across conditions?
fittedPars = array(NA, dim = c(5, nsubjs))
for (j in 1:nsubjs){
  print(j)
  #where does the conditions variable play its part?
  tmpdat = dat[dat$subjs == j & !is.na(dat$rt), ]
  outs = optim(par = startPar, logLikelihoodFunction, 
               parameterNames = parameterNames, 
               rt = tmpdat$rt, 
               resp = tmpdat$responses,
               conditions = tmpdat$conds)
    
    fittedPars[,j] = outs$par
  }

fittedPars
mns = apply(fittedPars, 1, mean)
mns
sds = apply(fittedPars, 1, sd)
sds

#Does drift rate vary across conditions?

#look at fit
nconditions = length(unique(dat$conds))

#have to choose what to look at - can start simple
obsChoiceProp = tapply(dat$responses == 2, list(dat$conds, dat$subjs), mean, na.rm =T)
obsUpMRT = tapply(dat$rt[dat$responses == 2], list(dat$conds[dat$responses == 2], dat$subjs[dat$responses == 2]), mean, na.rm = T)
obsLowMRT = tapply(dat$rt[dat$responses == 1], list(dat$conds[dat$responses == 1], dat$subjs[dat$responses == 1]), mean, na.rm = T)

matplot(x = 1:3, y = obsChoiceProp, pch = 1, col = 1, xlim = c(0.5,3.5), ylim = c(0,1))
matplot(x = 1:3, y = obsUpMRT, pch = 1, col = 1, xlim = c(0.5,3.5), ylim = c(0.3,2.5))
matplot(x = 1:3, y = obsLowMRT, pch = 1, col = 1, xlim = c(0.5,3.5), ylim = c(0.3,2.5))


#now get the same predicted data from the model
nsims = 1000
predUpMRT = predLowMRT = predChoiceProp = array(NA, dim = c(nconditions, nsubjs))

for (j in 1:nsubjs){
  tmpPars = fittedPars[,j]
  names(tmpPars) = parameterNames
  
  a = tmpPars[grep("a",parameterNames)]
  if (length(a) < nconditions) a = rep(a, nconditions)
  v = tmpPars[grep("v",parameterNames)]
  if (length(v) < nconditions) v = rep(v, nconditions)
  t0 = tmpPars[grep("t0",parameterNames)]
  if (length(t0) < nconditions) t0 = rep(t0, nconditions)
  
  #You need to put in the code to simulate from the model, and create the predicted data patterns
  
}

matplot(x = 1:3 - 0.1, y = obsChoiceProp, pch = 1, col = 1, xlim = c(0.5,3.5), ylim = c(0,1))
matpoints(x = 1:3 + 0.1, y = predChoiceProp, pch = 16,col = rgb(0,0,0,0.5))
matplot(obsChoiceProp - predChoiceProp, pch = 16, col = rgb(0,0,0,0.3))
abline(h = 0)

matplot(x = 1:3 - 0.1, y = obsUpMRT, pch = 1, col = 1, xlim = c(0.5,3.5), ylim = c(0.3,2.5))
matpoints(x = 1:3 + 0.1, y = predUpMRT, pch = 16,col = rgb(0,0,0,0.5))
matplot(obsUpMRT - predUpMRT, pch = 16, col = rgb(0,0,0,0.3))
abline(h = 0)

matplot(x = 1:3 - 0.1, y = obslowMRT, pch = 1, col = 1, xlim = c(0.5,3.5), ylim = c(0.3,2.5))
matpoints(x = 1:3 + 0.1, y = predLowMRT, pch = 16,col = rgb(0,0,0,0.5))
matplot(obsLowMRT - predLowMRT, pch = 16, col = rgb(0,0,0,0.3))
abline(h = 0)



### Try a different model ###

#change the code above so that a model with only boundary separation varies across conditions

#how does that do?

### Let's introduce start-point bias ###

logLikelihoodFunction=function(parameters,parameterNames,
                               rt,resp,conditions) {
  
  nconditions = length(unique(conditions))
  
  names(parameters) = parameterNames
  out=0
  
  a = parameters[grep("a",parameterNames)]
  if (length(a) < nconditions) a = rep(a, nconditions)
  v = parameters[grep("v",parameterNames)]
  if (length(v) < nconditions) v = rep(v, nconditions)
  t0 = parameters[grep("t0",parameterNames)]
  if (length(t0) < nconditions) t0 = rep(t0, nconditions)
  z = parameters[grep("z",parameterNames)]
  if (length(z) < nconditions) z = rep(z, nconditions)
  
  for (i in 1:nconditions){
    out = out + sum(log(pmax(
      ddiffusion(rt = rt[conditions == i], response = resp[conditions == i],
                 a = a[i], v = v[i], t0 = t0[i], z = z[i] * a[i], 
                 s = 0.1), 1e-10)))
  }
  -out
  
}

#going to have to fix this
startPar = c(0.05, 0.1, 0.15, 0.2, 0.15)
parameterNames = c("a1", "a2", "a3", "v", "t0")

fittedPars = array(NA, dim = c(5, nsubjs))
for (j in 1:nsubjs){
  print(j)
  
  tmpdat = dat[dat$subjs == j & !is.na(dat$rt), ]
  outs = optim(par = startPar, logLikelihoodFunction, 
               parameterNames = parameterNames, 
               rt = tmpdat$rt, 
               resp = tmpdat$responses,
               conditions = tmpdat$conds)
  
  fittedPars[,j] = outs$par
}

mns = apply(fittedPars, 1, mean)
sds = apply(fittedPars, 1, sd)

par(mar = c(4,4,1,1))
plot(mns, ylim = c(-1,1), axes = F, xlab = '', ylab = '')
axis(side = 1, at = 1:6, c('a1','a2','a3','v','t0','z')); axis(side = 2); box()
arrows(x0 = 1:6, y0 = mns-sds/sqrt(nsubjs), y1 = mns + sds/sqrt(nsubjs), code = 3, angle = 90, length = 0.05)

#check the predicted data patterns
for (j in 1:nsubjs){
  tmpPars = fittedPars[,j]
  names(tmpPars) = parameterNames
  
  a = tmpPars[grep("a",parameterNames)]
  if (length(a) < nconditions) a = rep(a, nconditions)
  v = tmpPars[grep("v",parameterNames)]
  if (length(v) < nconditions) v = rep(v, nconditions)
  t0 = tmpPars[grep("t0",parameterNames)]
  if (length(t0) < nconditions) t0 = rep(t0, nconditions)
  z = tmpPars[grep("z",parameterNames)]
  if (length(z) < nconditions) z = rep(z, nconditions)
  
  #simulate from the model and calculate the needed data patterns  
  
}

matplot(x = 1:3 - 0.1, y = predChoiceProp, pch = 1, col = 1, xlim = c(0.5,3.5), ylim = c(0,1))
matpoints(x = 1:3 + 0.1, y = obsChoiceProp, pch = 16,col = rgb(0,0,0,0.5))
matplot(obsChoiceProp - predChoiceProp, pch = 16, col = rgb(0,0,0,0.3))
abline(h = 0)

matplot(x = 1:3 - 0.1, y = predUpMRT, pch = 1, col = 1, xlim = c(0.5,3.5), ylim = c(0.3,2.5))
matpoints(x = 1:3 + 0.1, y = obsUpMRT, pch = 16,col = rgb(0,0,0,0.5))
matplot(obsUpMRT - predUpMRT, pch = 16, col = rgb(0,0,0,0.3))
abline(h = 0)

matplot(x = 1:3 - 0.1, y = predLowMRT, pch = 1, col = 1, xlim = c(0.5,3.5), ylim = c(0.3,2.5))
matpoints(x = 1:3 + 0.1, y = obsLowMRT, pch = 16,col = rgb(0,0,0,0.5))
matplot(obsLowMRT - predLowMRT, pch = 16, col = rgb(0,0,0,0.3))
abline(h = 0)


#################################
# can always check if drift changes if both boundary and drift vary


fittedPars = array(NA, dim = c(8, nsubjs))
for (j in 1:nsubjs){
  print(j)
  
  tmpdat = dat[dat$subjs == j & !is.na(dat$rt), ]
  outs = optim(par = startPar, logLikelihoodFunction, 
               parameterNames = parameterNames, 
               rt = tmpdat$rt, 
               resp = tmpdat$responses,
               conditions = tmpdat$conds)
  
  fittedPars[,j] = outs$par
}

mns = apply(fittedPars, 1, mean)
sds = apply(fittedPars, 1, sd)

par(mar = c(4,4,1,1))
plot(mns, ylim = c(-1,1), axes = F, xlab = '', ylab = '')
axis(side = 1, at = 1:8, c('a1','a2','a3','v1','v2','v3','t0','z')); axis(side = 2); box()
arrows(x0 = 1:8, y0 = mns-sds/sqrt(nsubjs), y1 = mns + sds/sqrt(nsubjs), code = 3, angle = 90, length = 0.05)


