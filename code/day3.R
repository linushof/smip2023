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

allV=seq(0,3,0.2)
allA=seq(0.2,3,0.2)

MCRT=array(NA,c(length(allV),length(allA)))
MERT=array(NA,c(length(allV),length(allA)))
ErrorProp=array(NA,c(length(allV),length(allA)))

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





