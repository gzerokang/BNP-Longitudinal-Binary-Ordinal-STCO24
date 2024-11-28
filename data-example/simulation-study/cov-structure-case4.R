library(mvtnorm)
library(MCMCpack)
library(BayesLogit)
library(CholWishart)
library(Matrix)
library(transport)
library(corrplot)
library(psych)
library(data.table)



expit=function(x){
  exp(x)/(1+exp(x))
}

#### variance covariance structure
## compound symmetry
var_cov_str=function(d,rho=0.4,sigma=1){
  sigma*ifelse(d==0,1,rho)
}



n=100 ## number of subjects
grid.T=seq(0,10,by=1/3) ## number of observation per subject

## distance matrix for prediction
D.mtx.all=dist(grid.T,diag=TRUE,upper=TRUE)
D.mtx.all=as.matrix(D.mtx.all)


## time index with observations 
T.idx=seq(1,length(grid.T),by=3) ## all time index with data
T.all=length(T.idx)

## distance matrix of the data
D.mtx=dist(grid.T[T.idx],diag=TRUE,upper=TRUE)
D.mtx=as.matrix(D.mtx)

## from some function
signal_fun=function(t){
  0.1+2*sin(0.5*t)+cos(0.5*t)
}

set.seed(1)
mean.vec.true=signal_fun(grid.T[T.idx]) ## 
sigma.epsilon.true=0.1 ## signal to noise ratio at 10
cov.mtx.true=var_cov_str(D.mtx)
noise=rmvt(n,df=5,sigma=cov.mtx.true) ## signal process is a TP
noise=noise+matrix(rnorm(n*length(grid.T[T.idx]),
                         mean=0,sd=sqrt(sigma.epsilon.true)),nrow=n) ## add a noise
Z.true=t(apply(noise, MARGIN=1, function(vec){vec+mean.vec.true})) 
true.prob.curve.mtx=expit(Z.true)

## no missing in this case
idx.mtx=matrix(1,nrow=n,ncol=T.all)

y.dat=matrix(0,nrow=n,ncol=T.all)

set.seed(1)
for (i in 1:n) {
  for (t in 1:T.all) {
    prob.cur=true.prob.curve.mtx[i,t]
    y.dat[i,t]=rbinom(1,1,prob=prob.cur)
  }
}


## visualize the true covariance structure

library(colorRamps)
library(leaflet)
library(fields)

myPalette=colorRampPalette(rev(RColorBrewer::brewer.pal(11,"Spectral")))
colours=myPalette(100)

cov.mtx.true.all=var_cov_str(D.mtx.all)+diag(sigma.epsilon.true,nrow=nrow(D.mtx.all)) 

image(grid.T,grid.T,cov.mtx.true.all,col=colours,xlab='Distance',
      ylab='Distance')


########################## Fit the model ################################
## loading functions
source("code/BinaryFunAll.R")

## setup
p=as.numeric(length(T.idx)) ## |\tau|
id.miss=seq(1:n)[rowSums(idx.mtx)!=p]

## prior hyperparameters
a.epsilon=0.4
b.epsilon=1
a.mu=0
b.mu=10
a.sigma=5
b.sigma=10
a.rho=6.0
b.rho=8.4
a.nu=4
b.nu=20

## MCMC iterations
nsim=30000
nburn=10000
nthin=4
nkeep=(nsim-nburn)/nthin

## matrix to store parameters
Z.mtx=matrix(0,nrow=p*n,ncol=nkeep)
sigma.epsilon.vec=rep(0,nkeep)
mu.mtx=matrix(0,nrow=p,ncol=nkeep)
Sigma.mtx=matrix(0,nrow=p*p,ncol=nkeep)
mu0.vec=rep(0,nkeep)
sigma.vec=rep(0,nkeep)
rho.vec=rep(1,nkeep)
nu.vec=rep(0,nkeep)

## initial values
latent.Z.all.cur=rep(1,p*n)
xi.cur=matrix(1,nrow=n,ncol=p)
Z.all.cur=rep(1,p*n)
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
                                                                  id.sub,idx.mtx,y.dat)
    
  }
  
  ## update latent \xi
  xi.cur=t(sapply(1:n, function(id.sub){sample_xi(latent.Z.all.cur,id.sub,idx.mtx)}))
     
  # ## update \tilde{Z}_i
  for (id.sub in 1:n) {
    if(id.sub %in% id.miss){
       Z.all.cur[((id.sub-1)*p+1):(id.sub*p)]=sample_Z_sub_ub(mu.cur,Sigma.cur,
                                                              sigma.epsilon.cur,latent.Z.all.cur,
                                                              Z.all.cur,id.sub,idx.mtx)
     }else{
       Z.all.cur[((id.sub-1)*p+1):(id.sub*p)]=sample_Z_sub_balance(mu.cur,Sigma.cur.inv,
                                                                   sigma.epsilon.cur,
                                                                   latent.Z.all.cur,id.sub,idx.mtx)
     }
  }
  
  diff.cur=cal_diff(latent.Z.all.cur,Z.all.cur,idx.mtx)
  # ## update \sigma_{\epsilon}
  sigma.epsilon.cur=sample_sigma_epsilon(a.epsilon,b.epsilon,diff.cur,idx.mtx)
   
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
  rho.cur=sample_rho(Sigma.cur,sigma.cur,a.rho,b.rho,nu.cur,length.rho=25)
  
  ## sample nu
  nu.cur=sample_nu(mu.cur,mu0.cur,Sigma.cur,sigma.cur,rho.cur,
                   a.nu,b.nu)
  kappa.cur=1/(nu.cur-3)
   
   
  ## store samples
  id.save=i-nburn
  if((id.save%/%nthin)>0 & (id.save%%nthin)==0){
    Z.mtx[,(id.save%/%nthin)]=Z.all.cur
    sigma.epsilon.vec[(id.save%/%nthin)]=sigma.epsilon.cur
    mu.mtx[,(id.save%/%nthin)]=as.vector(mu.cur)
    Sigma.mtx[,(id.save%/%nthin)]=as.vector(Sigma.cur)
    mu0.vec[(id.save%/%nthin)]=mu0.cur
    sigma.vec[(id.save%/%nthin)]=sigma.cur
    rho.vec[(id.save%/%nthin)]=rho.cur
    nu.vec[(id.save%/%nthin)]=nu.cur
  }
  
  
  
  ## print
  if(i %% 1000==0){
    print(i)
  }
  
}  
)




## check convergence

par(mfrow=c(2,3))
plot(sigma.epsilon.vec,type='l',xlab=expression(sigma[epsilon]^2),ylab='')
plot(mu0.vec,type='l',xlab=expression(mu[0]),ylab='')
plot(sigma.vec,type='l',xlab=expression(sigma[r]^2),ylab='')
plot(rho.vec,type='l',xlab=expression(rho),ylab='')
plot(nu.vec,type='l',xlab=expression(nu),ylab='')



## posterior inference of covariance structure for the latent process
## estimate the covariogram

post.varcov.sample.mtx=matrix(0,nrow=length(grid.T),ncol=nkeep)


for (i.s in 1:nkeep){
  set.seed(i.s)
  post.varcov.sample.mtx[,i.s]=post_cov_sample_signal(Dist.mtx=D.mtx.all,
                                                      sigma.sample=sigma.vec[i.s],
                                                      rho.sample=rho.vec[i.s])
}



mean.varcov.sample=rowMeans(post.varcov.sample.mtx)
upper.varcov.sample=apply(post.varcov.sample.mtx, MARGIN=1, function(vec){unname(quantile(vec,0.975))})
lower.varcov.sample=apply(post.varcov.sample.mtx, MARGIN=1, function(vec){unname(quantile(vec,0.025))})



## plot covariogram

plot(grid.T,var_cov_str(D.mtx.all)[1,],
     xlim=c(min(grid.T)-0.1,max(grid.T)+0.1),
     ylim=c(0,upper.varcov.sample[1]+0.1),
     xlab='Distance',ylab='Covariance',type='l',lty=1,col='#4DAF4A')
lines(grid.T,mean.varcov.sample,lty=2,col='#E41A1C')
polygon(c(grid.T,rev(grid.T)),
        c(upper.varcov.sample,rev(lower.varcov.sample)),
        col=adjustcolor("#A6CEE3", alpha.f = 0.25),border=NA)


