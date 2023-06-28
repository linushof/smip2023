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
