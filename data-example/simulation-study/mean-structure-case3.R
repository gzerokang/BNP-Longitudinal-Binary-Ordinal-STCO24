library(mvtnorm)
library(MCMCpack)
library(BayesLogit)
library(CholWishart)
library(Matrix)




expit=function(x){
  exp(x)/(1+exp(x))
}


power_exp=function(d,rho=1,sigma=1,power=2){
  sigma*exp(-(d/rho)^power)
}

n=30 ## number of subjects
grid.T=seq(0,30,by=1/3) ## number of observation per subject

## distance matrix for prediction
D.mtx.all=dist(grid.T,diag=TRUE,upper=TRUE)
D.mtx.all=as.matrix(D.mtx.all)


## generate unbalanced data
T.idx=seq(1,length(grid.T),by=3) ## all time index
T.all=length(T.idx)

## distance matrix of the data
D.mtx=dist(grid.T[T.idx],diag=TRUE,upper=TRUE)
D.mtx=as.matrix(D.mtx)

## from some function

signal_fun1=function(t){
  0.1+2*sin((1/4)*t)+cos(0.25*t)
}

signal_fun2=function(t){
  0.3+3*sin(0.5*t)+cos((1/3)*t)
}

set.seed(1)
mean.vec.true1=signal_fun1(grid.T[T.idx])
cov.mtx.true1=power_exp(D.mtx)
noise1=rmvt(n,df=5,sigma=1/3*cov.mtx.true1)
Z.true1=t(apply(noise1, MARGIN=1, function(vec){vec+mean.vec.true1})) ## continuous signal

mean.vec.true2=signal_fun2(grid.T[T.idx])
cov.mtx.true2=power_exp(D.mtx,power=1)
noise2=rmvnorm(n,rep(0,T.all),cov.mtx.true2)
Z.true2=t(apply(noise2, MARGIN=1, function(vec){vec+mean.vec.true2})) ## continuous signal

## adding the measurement error
set.seed(1)
measure.error1=matrix(rnorm(n*length(grid.T[T.idx]),mean=0,sd=sqrt(0.1)),nrow=n)
latent.process.mtx1=Z.true1+measure.error1
true.prob.curve.mtx1=expit(latent.process.mtx1)


## adding the measurement error
set.seed(123)
measure.error2=matrix(rnorm(n*length(grid.T[T.idx]),mean=0,sd=sqrt(0.1)),nrow=n)
latent.process.mtx2=Z.true2+measure.error2
true.prob.curve.mtx2=pnorm(latent.process.mtx2)

## missing pattern 
idx.mtx1=matrix(1,nrow=n,ncol=T.all)
miss.proportion1=0.75

## miss.proportion1=0.5 
## miss.proportion1=0.9 

for (i.sub in 1:n) {
  for (tau.time in 1:T.all) {
    set.seed(i.sub*tau.time)
    u=runif(1)
    if(u>miss.proportion1){
      idx.mtx1[i.sub,tau.time]=0
    }
  }
}

## check
## sum(idx.mtx1==0)/930

y.dat1=matrix(0,nrow=n,ncol=T.all)
mix.weight=0.5

set.seed(1)
for (i in 1:n) {
  for (t in 1:T.all) {
    if(idx.mtx1[i,t]==1){
      config=runif(1)
      if(config<=mix.weight){
        prob.cur1=true.prob.curve.mtx1[i,t]
      }else{
        prob.cur1=true.prob.curve.mtx2[i,t]
      }
      y.dat1[i,t]=rbinom(1,1,prob=prob.cur1)
    }else{
      y.dat1[i,t]=NA
    }
  }
}




############### Fit the model with the first data #########################
source("code/BinaryFunAll.R")


## setup
p=as.numeric(length(T.idx)) ## |\tau|
id.miss=seq(1:n)[rowSums(idx.mtx1)!=p]

## prior hyperparameters
a.epsilon=2
b.epsilon=0.1
a.mu=0
b.mu=10
a.sigma=4
b.sigma=4
a.rho=0.1
b.rho=2
a.nu=4
b.nu=20

## MCMC iterations
nsim=30000
nburn=10000
nthin=4
nkeep=(nsim-nburn)/nthin

## matrix to store parameters
Z.mtx1=matrix(0,nrow=p*n,ncol=nkeep)
sigma.epsilon.vec1=rep(0,nkeep)
mu.mtx1=matrix(0,nrow=p,ncol=nkeep)
Sigma.mtx1=matrix(0,nrow=p*p,ncol=nkeep)
mu0.vec1=rep(0,nkeep)
sigma.vec1=rep(0,nkeep)
rho.vec1=rep(1,nkeep)
nu.vec1=rep(0,nkeep)


## initial values
latent.Z.all.cur=rep(1,p*n)
xi.cur=matrix(1,nrow=n,ncol=p)
Z.all.cur=rep(1,p*n)
diff.cur=rep(0,n)
sigma.epsilon.cur=0.1
mu.cur=rep(0,p)
Sigma.cur=diag(p)
Sigma.cur.inv=diag(p)
mu0.cur=1
sigma.cur=1
rho.cur=1
nu.cur=4
kappa.cur=1/(nu.cur-3)



################# The loop ######################
system.time(
  for (i in 1:nsim) {
    set.seed(i)
    
    ## update latent \mathcal{Z}_i, \xi, and \tilde{Z}_i for each subject
    for (id.sub in 1:n) {
      ## update \mathcal{Z}_i
      latent.Z.all.cur[((id.sub-1)*p+1):(id.sub*p)]=sample_Z_latent(xi.cur,Z.all.cur,
                                                                    sigma.epsilon.cur,
                                                                    id.sub,idx.mtx1,y.dat1)
      
    }
    
    ## update latent \xi
    xi.cur=t(sapply(1:n, function(id.sub){sample_xi(latent.Z.all.cur,id.sub,idx.mtx1)}))
    
    # ## update \tilde{Z}_i
    for (id.sub in 1:n) {
      if(id.sub %in% id.miss){
        Z.all.cur[((id.sub-1)*p+1):(id.sub*p)]=sample_Z_sub_ub(mu.cur,Sigma.cur,
                                                               sigma.epsilon.cur,latent.Z.all.cur,
                                                               Z.all.cur,id.sub,idx.mtx1)
      }else{
        Z.all.cur[((id.sub-1)*p+1):(id.sub*p)]=sample_Z_sub_balance(mu.cur,Sigma.cur.inv,
                                                                    sigma.epsilon.cur,
                                                                    latent.Z.all.cur,id.sub,idx.mtx1)
      }
    }
    
    diff.cur=cal_diff(latent.Z.all.cur,Z.all.cur,idx.mtx1)
    # ## update \sigma_{\epsilon}
    sigma.epsilon.cur=sample_sigma_epsilon(a.epsilon,b.epsilon,diff.cur,idx.mtx1)
    
    # ## update mu
    mu.cur=sample_mu(mu0.cur,Sigma.cur,Z.all.cur,kappa.cur)
    mu.cur=matrix(mu.cur,ncol=1)
    
    ## update Sigma
    Sigma.cur=sample_Sigma(mu0.cur,Z.all.cur,rho.cur,sigma.cur,
                           nu.cur,kappa.cur)
    Sigma.cur.inv=chol2inv(chol(Sigma.cur))
    
    ## sample mu0
    mu0.cur=sample_mu0(mu.cur,Sigma.cur.inv,a.mu,b.mu)
    
    ## sample sigma
    sigma.cur=sample_sigma(rho.cur,Sigma.cur.inv,a.sigma,b.sigma,
                           nu.cur)
    
    ## sample rho
    rho.cur=sample_rho(Sigma.cur,sigma.cur,a.rho,b.rho,nu.cur,length.rho=20)
    
    ## sample nu
    nu.cur=sample_nu(mu.cur,mu0.cur,Sigma.cur,sigma.cur,rho.cur,
                     a.nu,b.nu)
    kappa.cur=1/(nu.cur-3)
    
    
    ## store samples
    id.save=i-nburn
    if((id.save%/%nthin)>0 & (id.save%%nthin)==0){
      Z.mtx1[,(id.save%/%nthin)]=Z.all.cur
      sigma.epsilon.vec1[(id.save%/%nthin)]=sigma.epsilon.cur
      mu.mtx1[,(id.save%/%nthin)]=as.vector(mu.cur)
      Sigma.mtx1[,(id.save%/%nthin)]=as.vector(Sigma.cur)
      mu0.vec1[(id.save%/%nthin)]=mu0.cur
      sigma.vec1[(id.save%/%nthin)]=sigma.cur
      rho.vec1[(id.save%/%nthin)]=rho.cur
      nu.vec1[(id.save%/%nthin)]=nu.cur
    }
    
    
    
    ## print
    if(i %% 1000==0){
      print(i)
    }
    
  }  
)

par(mfrow=c(2,3))
plot(sigma.epsilon.vec1,type='l',xlab=expression(sigma[epsilon]^2),ylab='')
plot(mu0.vec1,type='l',xlab=expression(mu[0]),ylab='')
plot(sigma.vec1,type='l',xlab=expression(sigma[r]^2),ylab='')
plot(rho.vec1,type='l',xlab=expression(rho),ylab='')
plot(nu.vec1,type='l',xlab=expression(nu),ylab='')

## posterior inference of probability response curve
## for in sample and out of sample

post.Z.sample1=matrix(0,nrow=n*length(grid.T),ncol=nkeep)
post.prob.sample1=matrix(0,nrow=n*length(grid.T),ncol=nkeep)
post.prob.sample.new=matrix(0,nrow=length(grid.T),ncol=nkeep)

system.time(
  for (i.s in 1:nkeep) {
    set.seed(i.s)
    ## inference the signal process
    post.sample.z=pred_Z_sample(rho.vec1[i.s],sigma.vec1[i.s],
                                mu0.vec1[i.s],Z.mtx1[,i.s],
                                nu.vec1[i.s])
    post.Z.sample1[,i.s]=post.sample.z
    ## inference the probability response curve
    post.prob.sample1[,i.s]=pred_prob_sample_delta_exact(post.sample.z,sigma.epsilon.vec1[i.s])
    post.prob.sample.new[,i.s]=pred_prob_sample_new(mu.mtx1[,i.s],Sigma.mtx1[,i.s],
                                              rho.vec1[i.s],sigma.vec1[i.s],
                                              mu0.vec1[i.s],nu.vec1[i.s],
                                              sigma.epsilon.vec1[i.s])
    if(i.s %% 1000==0){
      print(i.s)
    }
  }
)

## plot point and interval estimate of the probability response curve
mean.prob.all1=rowMeans(post.prob.sample1)
upper.prob.all1=apply(post.prob.sample1, MARGIN=1, function(vec){unname(quantile(vec,0.975))})
lower.prob.all1=apply(post.prob.sample1, MARGIN=1, function(vec){unname(quantile(vec,0.025))})


mean.prob.all.mtx1=matrix(mean.prob.all1,nrow=n,byrow=TRUE)
upper.prob.all.mtx1=matrix(upper.prob.all1,nrow=n,byrow=TRUE)
lower.prob.all.mtx1=matrix(lower.prob.all1,nrow=n,byrow=TRUE)


signal.true=mix.weight*signal_fun1(grid.T)+(1-mix.weight)*signal_fun2(grid.T)


## plot

par(mfrow=c(1,3),mar = c(4, 4, 1, 1))
set.seed(42)
sample.select=sample(1:n,size=3)
for (i in sample.select) {
  plot(grid.T[T.idx[idx.mtx1[i,]==1]],true.prob.curve.mtx1[i,idx.mtx1[i,]==1],xlim=c(min(grid.T)-0.1,max(grid.T)+0.1),
       ylim=c(-0.1,1.1),pch=4,
       xlab='T',ylab='Pr(Y=1)',type='p',col='gray48')
  for (j in 0:30) {
    abline(v=j,col='gray94',lty=3)
  }
  lines(grid.T,expit(signal.true),col='#4DAF4A',lty=1)
  points(grid.T[T.idx],y.dat1[i,],col='#FF7F00',pch=16)
  lines(grid.T,mean.prob.all.mtx1[i,],lty=2,col='#E41A1C')
  polygon(c(grid.T,rev(grid.T)),c(upper.prob.all.mtx1[i,],rev(lower.prob.all.mtx1[i,])),
          col=adjustcolor("#A6CEE3", alpha.f = 0.25),border=NA)
}


## RMSE of signal process
RMSE.signal.full.75=apply(post.Z.sample1, MARGIN=2, function(vec){
  mtx.cur=matrix(vec,nrow=n,byrow=TRUE)
  error.row=apply(mtx.cur, MARGIN=1, function(sign){mean((sign-signal.true)^2)})
  sqrt(mean(error.row))
})


## RMSE of probability
RMSE.prob.full.75=apply(post.prob.sample1, MARGIN=2, function(vec){
  mtx.cur=matrix(vec,nrow=n,byrow=TRUE)
  error.prob.vec=rep(0,n)
  for (id.sub in 1:n) {
    idx.with=idx.mtx1[id.sub,]==1
    T.idx.with=T.idx[idx.with]
    prob.sample.use=mtx.cur[id.sub,T.idx.with]
    true.prob.use=mix.weight*true.prob.curve.mtx1[id.sub,idx.with]+(1-mix.weight)*true.prob.curve.mtx2[id.sub,idx.with]
    error.prob.vec[id.sub]=mean((prob.sample.use-true.prob.use)^2)
  }
  sqrt(mean(error.prob.vec))
})


summary(RMSE.prob.full.75)



mean.prob.new=rowMeans(post.prob.sample.new)
upper.prob.new=apply(post.prob.sample.new, MARGIN=1, function(vec){unname(quantile(vec,0.975))})
lower.prob.new=apply(post.prob.sample.new, MARGIN=1, function(vec){unname(quantile(vec,0.025))})


plot(grid.T,expit(signal.true),xlim=c(min(grid.T)-0.1,max(grid.T)+0.1),
     ylim=c(-0.1,1.1),pch=4,
     xlab='T',ylab='Pr(Y=1)',type='l',col='#4DAF4A',lwd=3)
for (id.sub in 1:n) {
  lines(grid.T,mean.prob.all.mtx1[id.sub,],lty=1,col='gray64')
}
lines(grid.T,mean.prob.new,lty=2,col='#E41A1C',lwd=3)
polygon(c(grid.T,rev(grid.T)),c(upper.prob.new,rev(lower.prob.new)),
        col=adjustcolor("#A6CEE3", alpha.f = 0.25),border=NA)


