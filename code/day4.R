simDiffusion=function(N,v,as,ell,ter,z,s,dt,maxiter) {
  resp = rep(-1,N)
  rt = rep(-1,N)
  x = array(NA,c(N,maxiter+1))
  av = as * exp(-ell*(1:maxiter)) 
  x[,1] = z
  iter = 0
  keep = rep(TRUE,N)
  a = av[iter+1]
  while (iter < maxiter) {
    iter = iter + 1
    x[keep,iter+1] = x[keep,iter] + v*dt + s*sqrt(dt)*rnorm(sum(keep))
    if (any(x[keep,iter+1] > a)) {
      resp[keep & x[,iter+1]> a] = 2
      rt[keep & x[,iter+1]> a] = ter + iter*dt - dt/2
      keep[keep & x[,iter+1]> a] = FALSE
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
v=1
as=2
ell = .0002
ter=0.3
z=as/2
s=1
dt=0.001
maxiter=10000


begin1=proc.time()
tmp1=simDiffusion(N=N,v=v,as=as, ell = ell, ter=ter,z=z,s=s,dt=dt,maxiter=maxiter)
end1=proc.time()
print(end1-begin1)
cat("Correct Proportion: ",mean(tmp1$resp==2),"\n")
cat("Error Proportion: ",mean(tmp1$resp==1),"\n")
cat("Mean Correct RT: ",mean(tmp1$rt[tmp1$resp==2]),"\n")
cat("Mean Error RT: ",mean(tmp1$rt[tmp1$resp==1]),"\n")
hist(tmp1$rt[tmp1$resp==2])
plot(density(tmp1$rt[tmp1$resp==2]))
hist(tmp1$rt[tmp1$resp==1])
plot(density(tmp1$rt[tmp1$resp==1]))


begin2=proc.time()
tmp2=rdiffusion(n=N,a=a,v=v,t0=ter)
end2=proc.time()
print(end2-begin2)
tmp3=list(rt=tmp2$rt,resp=as.numeric(tmp2$response))
cat("Correct Proportion: ",mean(tmp3$resp==2),"\n")
cat("Error Proportion: ",mean(tmp3$resp==1),"\n")
cat("Mean Correct RT: ",mean(tmp3$rt[tmp3$resp==2]),"\n")
cat("Mean Error RT: ",mean(tmp3$rt[tmp3$resp==1]),"\n")
hist(tmp3$rt[tmp3$resp==2])
plot(density(tmp3$rt[tmp3$resp==2]))
hist(tmp3$rt[tmp3$resp==1])
plot(density(tmp3$rt[tmp3$resp==1]))