pacman::p_load("rtdists", "tidyverse")

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

