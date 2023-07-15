
#======== Causal Inference in CRT with interference ========#

library(lme4)
library(mvtnorm)
library(geepack)
#-------- parameters --------

# set.seed(123456)


#-------- simulation function ------
f.delt  = function(mu){
  delt0. = mu[2]-mu[1]
  delt1. = mu[4]-mu[3]
  delt.0 = mu[3]-mu[1]
  delt.1 = mu[4]-mu[2]
  c(delt0. , delt1. , delt.0, delt.1)
}

f.data = function(dat){
  M = length(unique(dat$clust))
  N = nrow(dat)
  ni = table(dat$clust)
  #-------- esitmate ------
  mu.GC = mu.IPW1 = my.IPW2 = rep(NA,4)
  
  #---- GC
  
  ml.y = formula(paste("y~",paste(paste("x",0:15,sep=""),collapse="+"),"+tr*s+ x14:tr +(1|clust)",sep =""))
  
  fit.y = lmer(ml.y, data = dat)  
  
  f.GC = function(trs, dat, fit.y){
    dat$tr= trs[1]
    dat$s = trs[2]
    dat$s.nj = trs[2]
    predict(fit.y, dat, re.form=~0)
  }
  tr.s = cbind( rep(0:1, c(2,2)),rep(0:1,2))
  
  mu.GC = colMeans(apply(tr.s, 1, f.GC, dat=dat, fit.y=fit.y))
  
  
  #---- IPW
  
  tr.clust = dat$tr[!duplicated(dat$clust)] 
  pred.pi = mean(tr.clust)
  
  ml.s = formula(paste("s~",paste(paste("x",0:15,sep=""),collapse="+"),"+tr+(1|clust)",sep =""))
  
  fit.s = glmer(ml.s, data = dat, family = binomial(link="probit"),
                control = glmerControl(check.conv.grad = .makeCC("warning", tol = 2e-1, relTol = NULL)))
  
  tau.s = predict(fit.s,dat, re.form=~0, type = "link")
  sig2.a = diag(summary(fit.s)$varcor$clust)  
  
  pred.s = rep(NA,M)
  s.clust = matrix(NA,M,2)
  y.clust = rep(NA,M)
  for (i in 1:M){
    indx.i = (dat$clust==i)
    dat.i = dat[indx.i,]
    s.clust[i,1] = all(dat.i$s==0)
    s.clust[i,2] = all(dat.i$s==1)
    
    sign.s = (2*dat.i$s-1)
    tau = tau.s[indx.i]*sign.s
    n = ni[i]
    sig.tau = matrix(sig2.a, n, n)*tcrossprod(sign.s,sign.s)
    diag(sig.tau) = 1+sig2.a
    
    pred.s[i] = pmvnorm(lower=-Inf, upper=tau, mean=rep(0, length(tau)), sigma=sig.tau)
    y.clust[i] = sum(dat.i$y)
    
  }
  
  rpr.pi = ifelse(tr.clust==1, pred.pi, 1-pred.pi)
  rpr.s = pred.s
  
  rpr = 1/rpr.pi/rpr.s
  y.ipw = y.clust*rpr
  
  f.ipw1.t1 = function(trs){
    tr = trs[1]
    s = trs[2]
    stemp = s+1
    ifelse(any((tr.clust==tr)&(s.clust[,stemp])), sum(y.ipw[(tr.clust==tr)&(s.clust[,stemp])])/N, NA)
  }
  mu.IPW1 = apply(tr.s, 1, f.ipw1.t1)
  
  fit.gee = geeglm(y~tr*s, data= dat, id=clust, family=gaussian("identity"),
                   weights = rep(rpr,ni),
                   corstr = "exchangeable") 
  # print(fit.gee)
  
  x.gee = data.frame(tr=tr.s[,1], s=tr.s[,2])
  
  mu.IPW2 = predict(fit.gee,x.gee)

  mu = cbind(mu.GC, mu.IPW1, mu.IPW2)
  return(mu)
}

dat=read.csv("DRINK_rate.csv")[,-1]
mu = f.data(dat)
delt = apply(mu,2,f.delt)

print(mu)
print(delt)

#---- bootstrap ----

f.boots.pall = function(func, dat){
  library(lme4)
  library(mvtnorm)
  library(geepack)
  
  M = length(unique(dat$clust))
  ind = sample(1:M, M, replace=T)
  
  dat.boots = c()
  for(i in 1:M){
    dat.i = dat[dat$clust==ind[i],]
    dat.i$clust = i
    dat.boots = rbind(dat.boots, dat.i)
  }
  

  mu.boots.try = try(func(dat.boots), silent = TRUE)
  if (inherits(mu.boots.try, "try-error")) { mu.boots=rep(NA,12)}else{
    mu.boots = mu.boots.try
  }
  return(c(mu.boots))
}

nboots = 500

star.sim = 1
end.sim = star.sim + nboots-1


#---bootstrapping

res=matrix(NA,nboots,12)
for(iboots in 1:nboots){
  res.mu[i,] =  f.boots.pall(func = f.data, dat)
  
}

f.each = function(x){
  y = matrix(x,ncol=3)
  f.all = function(mu){
    delt0. = mu[2]-mu[1]
    delt1. = mu[4]-mu[3]
    delt.0 = mu[3]-mu[1]
    delt.1 = mu[4]-mu[2]
    tau = mu[4]+mu[1]-mu[2]-mu[3]
    delt.. = mu[4]-mu[1]
    c(mu,delt0. , delt1. , delt.0, delt.1, tau, delt..)
  }
  dif.mu = c(apply(y,2,f.all))
}

#"Mean"
f.each(c(mu))
#"SD"
dif = t(apply(res.mu,1,f.each))
apply(dif,2,sd)
#"p-value"
2*(1-pnorm(abs(fmu[1,])/fmu[2,]))














  