cov_fun=function(d,rho=1,sigma=1){
   sigma*(1+sqrt(5)*d/rho+(5*d^2)/(3*rho^2))*exp(-sqrt(5)*d/rho)
}



data_processing=function(y.dat,C=2){
  y.dat.list=list()
  idx.mtx.list=list()
  for(j in 1:(C-1)){
    y.dat.list[[j]]=matrix(0,nrow=n,ncol=T.all)
    idx.mtx.list[[j]]=matrix(0,nrow=n,ncol=T.all)
  }
  
  for (i in 1:n) {
    for (t in 1:T.all) {
      if(is.na(y.dat[i,t])){
        for (j in 1:(C-1)) {
          y.dat.list[[j]][i,t]=NA
          idx.mtx.list[[j]][i,t]=0
        }
      }else{
        for (j in 1:(C-1)) {
          if(y.dat[i,t]>=j){
            idx.mtx.list[[j]][i,t]=1
            if(y.dat[i,t]==j){
              y.dat.list[[j]][i,t]=1
            }
          }
        }
      }
    }
  }
  
  return(c(y.dat.list,idx.mtx.list))
}

pool_grid=function(idx.mtx){
  idx.all=seq(1,ncol(idx.mtx),by=1)
  sort(unique(unlist(apply(idx.mtx,MARGIN = 1,function(vec){idx.all[vec==1]}))))
}


####################### MCMC functions #################################

sample_Z_latent=function(xi.cur,Z.all.cur,sigma.epsilon.cur,id.sub,idx.mtx,y.dat){
  ## prepare
  latent.Z.sub.sample=rep(0,p)
  idx.vec=seq(1,p,by=1)
  idx.sub=idx.vec[idx.mtx[id.sub,]==1] ## index with data
  Z.sub.cur=matrix(Z.all.cur,nrow=n,byrow=TRUE)[id.sub,]
  
  diag.vec=xi.cur[id.sub,idx.sub]+1/sigma.epsilon.cur
  cal.V.sub=diag(1/diag.vec,nrow=length(diag.vec))
  #cal.V.sub=chol2inv(chol(Omega.mtx+diag(1/sigma.epsilon.cur,length(idx.sub))))
  
  m.sub.part2=(1/sigma.epsilon.cur)*matrix(Z.sub.cur[idx.sub],ncol=1)
  lambda.vec=matrix(y.dat[id.sub,idx.sub]-1/2,ncol=1)
  m.sub=cal.V.sub%*%(lambda.vec+m.sub.part2)
  
  latent.Z.sub.sample.new=rmvnorm(1,m.sub,cal.V.sub)
  latent.Z.sub.sample[idx.sub]=latent.Z.sub.sample.new
  
  return(latent.Z.sub.sample)
}



sample_xi=function(latent.Z.all.cur,id.sub,idx.mtx){
  xi.sub.sample=rep(0,p)
  idx.vec=seq(1,p,by=1)
  idx.sub=idx.vec[idx.mtx[id.sub,]==1]
  latent.Z.sub.cur=matrix(latent.Z.all.cur,nrow=n,byrow=TRUE)[id.sub,idx.sub]
  xi.sample.new=sapply(1:length(latent.Z.sub.cur),function(idx){rpg.devroye(1,h=1,
                                                                     z=latent.Z.sub.cur[idx])})
  xi.sub.sample[idx.sub]=xi.sample.new
  return(xi.sub.sample)
}


cal_diff=function(latent.Z.all.cur,Z.all.cur,idx.mtx){
  latent.Z.all.cur.mtx=matrix(latent.Z.all.cur,nrow=n,byrow=TRUE)
  Z.all.cur.mtx=matrix(Z.all.cur,nrow=n,byrow=TRUE)
  diff.vec=sapply(1:n, function(idx){
    idx.select.vec=(idx.mtx[idx,]==1) ## index with data
    sum((latent.Z.all.cur.mtx[idx,idx.select.vec]-Z.all.cur.mtx[idx,idx.select.vec])^2)
  })
  sum(diff.vec)
}

sample_sigma_epsilon=function(a.epsilon,b.epsilon,diff.cur,idx.mtx){
  shape.sigma.epsilon=a.epsilon+sum(idx.mtx)/2
  rate.sigma.epsilon=b.epsilon+diff.cur/2
  sigma.epsilon.sample=1/rgamma(1,shape=shape.sigma.epsilon,
                                rate=rate.sigma.epsilon)
  return(sigma.epsilon.sample)
}

sample_Z_sub_balance=function(mu.cur,Sigma.cur.inv,sigma.epsilon.cur,latent.Z.all.cur,id.sub,idx.mtx){
  
  idx.vec=seq(1,p,by=1)
  idx.sub=idx.vec[idx.mtx[id.sub,]==1]
  
  I.mtx=diag((1/sigma.epsilon.cur),nrow=length(idx.sub))
  V.tilde.sub=chol2inv(chol(I.mtx+Sigma.cur.inv))
  
  mu.tilde2=Sigma.cur.inv%*%mu.cur
  latent.Z.sub.cur=matrix(latent.Z.all.cur,nrow=n,byrow=TRUE)[id.sub,]
  mu.tilde1=(1/sigma.epsilon.cur)*matrix(latent.Z.sub.cur[idx.sub],ncol=1)
  mu.tilde=V.tilde.sub%*%(mu.tilde1+mu.tilde2)
  
  Z.sub.sample=rmvnorm(1,mu.tilde,V.tilde.sub)
  return(Z.sub.sample)
  
}

sample_Z_sub_ub=function(mu.cur,Sigma.cur,sigma.epsilon.cur,latent.Z.all.cur,Z.all.cur,id.sub,idx.mtx){
  ## prepare
  Z.sub.sample=rep(0,length(T.idx))
  idx.vec=seq(1,length(T.idx),by=1)
  idx.sub=idx.vec[idx.mtx[id.sub,]==1]
  Sigma.tau.tau.inv=chol2inv(chol(Sigma.cur[idx.sub,idx.sub]))
  Bi.sub=Sigma.cur[-idx.sub,idx.sub]%*%Sigma.tau.tau.inv
  ui.sub=Bi.sub%*%matrix(mu.cur[idx.sub],ncol=1)-matrix(mu.cur[-idx.sub],ncol=1)
  Z.sub.cur=matrix(Z.all.cur,nrow=n,byrow=TRUE)[id.sub,]
  latent.Z.sub.cur=matrix(latent.Z.all.cur,nrow=n,byrow=TRUE)[id.sub,]
  
  ## update Z_i^*
  mean.sub.star=Bi.sub%*%matrix(Z.sub.cur[idx.sub],ncol=1)-ui.sub
  cov.mtx.sub.star=Sigma.cur[-idx.sub,-idx.sub]-Bi.sub%*%Sigma.cur[idx.sub,-idx.sub]
  cov.mtx.sub.star=forceSymmetric(cov.mtx.sub.star)
  cov.mtx.sub.star=as.matrix(cov.mtx.sub.star)
  Z.sub.star.sample=rmvnorm(1,mean.sub.star,cov.mtx.sub.star)
  Z.sub.sample[-idx.sub]=Z.sub.star.sample
  
  ## update Z_i
  I.mtx=diag(1/sigma.epsilon.cur,nrow=length(idx.sub))
  cov.mtx.sub.star.inv=chol2inv(chol(cov.mtx.sub.star))
  V.tilde=chol2inv(chol(I.mtx+Sigma.tau.tau.inv+
                          t(Bi.sub)%*%cov.mtx.sub.star.inv%*%Bi.sub))
  mu.tilde1=(1/sigma.epsilon.cur)*matrix(latent.Z.sub.cur[idx.sub],ncol=1)
  mu.tilde2=Sigma.tau.tau.inv%*%matrix(mu.cur[idx.sub],ncol=1)
  mu.tilde3=t(Bi.sub)%*%cov.mtx.sub.star.inv%*%(ui.sub+matrix(Z.sub.star.sample,ncol=1))
  mu.tilde=V.tilde%*%(mu.tilde1+mu.tilde2+mu.tilde3)
  Z.sub.in.sample=rmvnorm(1,mu.tilde,V.tilde)
  Z.sub.sample[idx.sub]=Z.sub.in.sample
  
  return(Z.sub.sample)
}



sample_mu=function(mu0.cur,Sigma.cur,Z.all.cur,kappa.cur){
  Z.all.mtx=matrix(Z.all.cur,ncol=n,byrow=FALSE)
  Z.all.mean=rowMeans(Z.all.mtx)
  mu0.vec=mu0.cur*matrix(1,nrow=p,ncol=1)
  mu.star=kappa.cur/(kappa.cur+n)*mu0.vec+n/(kappa.cur+n)*Z.all.mean
  kappa.star=n+kappa.cur
  
  mu.sample=rmvnorm(1,mean=mu.star,sigma=Sigma.cur/kappa.star)
  return(mu.sample)
}

sample_Sigma=function(mu0.cur,Z.all.cur,rho.cur,sigma.cur,nu.cur,kappa.cur){
  Z.all.mtx=matrix(Z.all.cur,ncol=n,byrow=FALSE)
  Z.all.mean=rowMeans(Z.all.mtx)
  Z.error=apply(Z.all.mtx-Z.all.mean,MARGIN = 2,function(x) x%*%t(x))
  S=apply(Z.error, MARGIN =1 , sum)
  S=matrix(S,nrow=p)
  nu.star=nu.cur+n
  df.param=nu.star+p-1
  mu0.vec=mu0.cur*matrix(1,nrow=p,ncol=1)
  
  Psi.use=cov_fun(D.mtx,rho=rho.cur,sigma=sigma.cur)
  
  Psi.star=Psi.use+S+n*kappa.cur*(Z.all.mean-mu0.vec)%*%t(Z.all.mean-mu0.vec)/(kappa.cur+n)
  
  cov.mtx.sample=rInvWishart(n=1,df=df.param, Sigma=Psi.star)[,,1]
  return(cov.mtx.sample)
}

sample_mu0=function(mu.cur,Sigma.cur.inv,nu.cur,a.mu,b.mu){
  one.vec=matrix(1,nrow=p,ncol=1)
  factor=1/(nu.cur-3)
  b.mu.star=1/(factor*t(one.vec)%*%Sigma.cur.inv%*%one.vec+1/b.mu)
  a.mu.star=b.mu.star*(factor*t(one.vec)%*%Sigma.cur.inv%*%mu.cur+a.mu/b.mu)
  mu0.sample=rnorm(1,mean=a.mu.star,sd=sqrt(b.mu.star))
  return(mu0.sample)
}

sample_sigma=function(rho.cur,Sigma.cur.inv,a.sigma,b.sigma,nu.cur){
  shape.sigma=a.sigma+(nu.cur+p-1)*p/2
  Psi.use=cov_fun(D.mtx,rho=rho.cur,sigma=1)
  rate.sigma=b.sigma+0.5*sum(diag(Psi.use%*%Sigma.cur.inv))
  sigma.sample=rgamma(1,shape=shape.sigma,rate=rate.sigma)
  return(sigma.sample)
}


sample_rho=function(Sigma.cur,sigma.cur,a.rho,b.rho,nu.cur,length.rho=30){
  rho.possible=seq(a.rho,b.rho,length.out=length.rho)
  log.prob.rho=sapply(1:length(rho.possible), function(idx){
    Psi.cur=cov_fun(D.mtx,rho=rho.possible[idx],sigma=sigma.cur)
    dInvWishart(Sigma.cur,nu.cur+p-1,Psi.cur)
  })
  prob.rho=exp(ifelse(log.prob.rho>=700,700,log.prob.rho))
  rho.idx=sample(1:length(rho.possible),1,replace=FALSE,prob=prob.rho/(sum(prob.rho)))
  
  rho.sample=rho.possible[rho.idx]
  return(rho.sample)
}


sample_nu=function(mu.cur,mu0.cur,Sigma.cur,sigma.cur,rho.cur,a.nu,b.nu,by.nu=1){
  nu.possible=seq(a.nu,b.nu,by=by.nu)
  Psi.cur=cov_fun(D.mtx,rho=rho.cur,sigma=sigma.cur)
  log.prob.nu=sapply(1:length(nu.possible), function(idx){
    mu.mean.vec.cur=rep(mu0.cur,length(mu.cur))
    mu.cov.mtx.cur=(nu.possible[idx]-3)*Sigma.cur
    log.lik.mu=dmvnorm(as.vector(mu.cur),mean=mu.mean.vec.cur,sigma=mu.cov.mtx.cur,log=TRUE)
    log.lik.Sigma=dInvWishart(Sigma.cur,nu.possible[idx]+p-1,Psi.cur)
    log.lik.mu+log.lik.Sigma
  })
  prob.nu=exp(ifelse(log.prob.nu>=700,700,log.prob.nu))
  nu.idx=sample(1:length(nu.possible),1,replace=FALSE,prob=prob.nu/(sum(prob.nu)))
  
  nu.sample=nu.possible[nu.idx]
  return(nu.sample)
}


############## posterior inference for covariance structure ###########


post_cov_sample_latent=function(Dist.mtx,sigma.epsilon.sample,sigma.sample,
                         rho.sample){
  cov.mtx.sample=cov_fun(Dist.mtx,rho=rho.sample,sigma=sigma.sample)
  #cov.mtx.sample=rInvWishart(n=1,df=nu.sample+nrow(Dist.mtx)-1, Sigma=cov.mtx)[,,1]
  var.mtx.sample=diag(sigma.epsilon.sample,nrow=nrow(Dist.mtx))
  varcov.mtx.sample=var.mtx.sample+cov.mtx.sample
  varcov.mtx.sample=varcov.mtx.sample[1,]
  return(varcov.mtx.sample)
}

post_cov_sample_signal=function(Dist.mtx,sigma.sample,rho.sample){
  cov.mtx.sample=cov_fun(Dist.mtx,rho=rho.sample,sigma=sigma.sample)
  cov.mtx.sample=cov.mtx.sample[1,]
  return(cov.mtx.sample)
}

## calculate Wasserstein distance
cal_wassd=function(Dist.mtx,nu.sample,sigma.sample,rho.sample,mu0.sample=0,nsim=500){
  cov.mtx.sample=cov_fun(Dist.mtx,rho=rho.sample,sigma=sigma.sample)
  sample.signal=rmvt(nsim,df=nu.sample,sigma=cov.mtx.sample)
  sample.signal.pp=pp(sample.signal)
  wasserstein(true.signal.pp,sample.signal.pp,p=2)
}



## posterior predictive for grid \boldsymbol{\tau}
pred_prob_sample_new_in=function(mu.sample,Sigma.sample,sigma.epsilon.sample){
  
  Cov.mtx.use=matrix(Sigma.sample,nrow=p,ncol=p)
  mean.vec.use=mu.sample
  Z.sample.in.new=rmvnorm(1,mean.vec.use,Cov.mtx.use)
  
  prob.sample.in.new=expit(Z.sample.in.new+rnorm(length(Z.sample.in.new),mean=0,sd=sqrt(sigma.epsilon.sample)))
  
  return(prob.sample.in.new)
}


## covariance
cal_COV_data=function(loc,dist,dat.mtx){
  if(dist==0){
    sub.vec=dat.mtx[,loc]
    true.cov.part=sub.vec-mean(sub.vec)
    true.cov=mean(true.cov.part*true.cov.part)
  }else{
    sub.mtx=dat.mtx[,c(loc,(loc+dist))]
    true.cov.part1=sub.mtx[,1]-mean(sub.mtx[,1])
    true.cov.part2=sub.mtx[,2]-mean(sub.mtx[,2])
    true.cov=mean(true.cov.part1*true.cov.part2)
  }
  return(true.cov)
}

cal_COV_data_vec=function(dist,dat.mtx){
  possible.loc=seq(1,(p-dist),by=1)
  COV.data.vec=sapply(1:length(possible.loc), function(idx){
    cal_COV_data(loc=possible.loc[idx],dist=dist,
                 dat.mtx=dat.mtx)
  })
  return(COV.data.vec)
}


## phi correlation
cal_phi_data=function(loc,dist,dat.mtx){
  sub.mtx=dat.mtx[,c(loc,(loc+dist))]
  p00=0
  p01=0
  p10=0
  p11=0
  for (row in 1:nrow(sub.mtx)) {
    bin.data.cur=sub.mtx[row,]
    if(bin.data.cur[1]==0 & bin.data.cur[2]==0){
      p00=p00+1
    }else if(bin.data.cur[1]==0 & bin.data.cur[2]==1){
      p01=p01+1
    }else if(bin.data.cur[1]==1 & bin.data.cur[2]==0){
      p10=p10+1
    }else{
      p11=p11+1
    }
  }
  phi.loc.dist=(p00*p11-p01*p10)/sqrt((p00+p10)*(p00+p01)*(p11+p01)*(p11+p10))
  return(phi.loc.dist)
}

cal_phi_data_vec=function(dist,dat.mtx){
  possible.loc=seq(1,(p-dist),by=1)
  phi.data.vec=sapply(1:length(possible.loc), function(idx){
    cal_phi_data(loc=possible.loc[idx],dist=dist,
                 dat.mtx=dat.mtx)
  })
  return(phi.data.vec)
}


## tet correlation
cal_tet_data=function(loc,dist,dat.mtx){
  sub.mtx=dat.mtx[,c(loc,(loc+dist))]
  p00=0
  p01=0
  p10=0
  p11=0
  for (row in 1:nrow(sub.mtx)) {
    bin.data.cur=sub.mtx[row,]
    if(bin.data.cur[1]==0 & bin.data.cur[2]==0){
      p00=p00+1
    }else if(bin.data.cur[1]==0 & bin.data.cur[2]==1){
      p01=p01+1
    }else if(bin.data.cur[1]==1 & bin.data.cur[2]==0){
      p10=p10+1
    }else{
      p11=p11+1
    }
  }
  tet.loc.dist=tetrachoric(matrix(c(p00,p01,p10,p11),2,2))$rho
  #tet.loc.dist.in=sqrt(p01*p10)/(sqrt(p00*p11)+sqrt(p01*p10))
  #tet.loc.dist=cospi(tet.loc.dist.in)
  return(tet.loc.dist)
}

cal_tet_data_vec=function(dist,dat.mtx){
  possible.loc=seq(1,(p-dist),by=1)
  tet.data.vec=sapply(1:length(possible.loc), function(idx){
    cal_tet_data(loc=possible.loc[idx],dist=dist,
                 dat.mtx=dat.mtx)
  })
  return(tet.data.vec)
}




################# prediction for subjects in sample ################

pred_Z_sample=function(rho.sample,sigma.sample,mu0.sample,Z.all.sample,nu.cur){
  Psi.mtx.post=cov_fun(D.mtx.all,rho=rho.sample,sigma=sigma.sample)
  mu0.vec.post=matrix(mu0.sample,nrow=length(grid.T),ncol=1)
  
  Z.mtx.sample=matrix(Z.all.sample,nrow=n,byrow=TRUE) ## each row is for each subject
  
  length.pred=length(grid.T)
  Z.out.sample=rep(0,n*length.pred)
  for (i.sub in 1:n) {
    Z.sample.sub=rep(0,length.pred)
    idx.sub=T.idx
    
    Psi.ii=Psi.mtx.post[idx.sub,idx.sub]
    Psi.ii.inv=chol2inv(chol(Psi.ii))
    Psi.oo=Psi.mtx.post[-idx.sub,-idx.sub]
    Psi.oi=Psi.mtx.post[-idx.sub,idx.sub]
    Psi.io=Psi.mtx.post[idx.sub,-idx.sub]
    
    tilde.Psi.oo=Psi.oo-Psi.oi%*%Psi.ii.inv%*%Psi.io
    #tilde.Psi.oo=forceSymmetric(tilde.Psi.oo)
    #tilde.Psi.oo=as.matrix(tilde.Psi.oo)
    
    mu.o=mu0.vec.post[-idx.sub]
    mu.i=mu0.vec.post[idx.sub]
    
    diff.i=matrix(Z.mtx.sample[i.sub,]-mu.i,ncol=1)
    
    tilde.mu.o=mu.o+Psi.oi%*%Psi.ii.inv%*%diff.i
    
    df.o=nu.cur+length(T.idx)
    
    S1=as.numeric(t(diff.i)%*%Psi.ii.inv%*%diff.i)
    fact.cov=(nu.cur+S1-2)/df.o ## (nu+S1-2)/(df.o-2)*(df.o-2)/df.o
    ## factor for scale matrix * factor for prediction
    
    pred.scale.mtx=fact.cov*tilde.Psi.oo
    pred.scale.mtx=(pred.scale.mtx+t(pred.scale.mtx))/2
    
    Z.out.sample.new=rmvt(1,sigma=pred.scale.mtx,df=df.o,delta=tilde.mu.o,
                          type='shifted')
    Z.sample.sub[-idx.sub]=Z.out.sample.new
    Z.sample.sub[idx.sub]=Z.mtx.sample[i.sub,]
    
    Z.out.sample[((i.sub-1)*length.pred+1):(i.sub*length.pred)]=Z.sample.sub
  }
  return(Z.out.sample)
}

expit_2nd_derv=function(x){
  2*(expit(x))^3-3*(expit(x))^2+expit(x)
}

pred_prob_sample_delta_approx=function(pred.z.sample,sigma.epsilon.sample){
  expit(pred.z.sample)+sigma.epsilon.sample*expit_2nd_derv(pred.z.sample)/2
}


pred_prob_sample_delta_exact=function(pred.z.sample,sigma.epsilon.sample){
  expit(pred.z.sample+rnorm(length(pred.z.sample),mean=0,sd=sqrt(sigma.epsilon.sample)))
}





pred_prob_sample_new=function(mu.sample,Sigma.sample,rho.sample,sigma.sample,mu0.sample,
                              nu.cur,sigma.epsilon.sample){
  Psi.mtx.post=cov_fun(D.mtx.all,rho=rho.sample,sigma=sigma.sample)
  mu0.vec.post=matrix(mu0.sample,nrow=length(grid.T),ncol=1)
  
  Cov.mtx.use=matrix(Sigma.sample,nrow=p,ncol=p)
  mean.vec.use=mu.sample
  Z.sample.in.new=rmvnorm(1,mean.vec.use,Cov.mtx.use)
  
  length.pred=length(grid.T)
  idx.sub=T.idx
  Z.sample.new=rep(0,length.pred)
  
  Psi.ii=Psi.mtx.post[idx.sub,idx.sub]
  Psi.ii.inv=chol2inv(chol(Psi.ii))
  Psi.oo=Psi.mtx.post[-idx.sub,-idx.sub]
  Psi.oi=Psi.mtx.post[-idx.sub,idx.sub]
  Psi.io=Psi.mtx.post[idx.sub,-idx.sub]
  
  tilde.Psi.oo=Psi.oo-Psi.oi%*%Psi.ii.inv%*%Psi.io
  #tilde.Psi.oo=forceSymmetric(tilde.Psi.oo)
  #tilde.Psi.oo=as.matrix(tilde.Psi.oo)
  
  mu.o=mu0.vec.post[-idx.sub]
  mu.i=mu0.vec.post[idx.sub]
  
  diff.i=matrix(Z.sample.in.new-mu.i,ncol=1)
  
  tilde.mu.o=mu.o+Psi.oi%*%Psi.ii.inv%*%diff.i
  
  df.o=nu.cur+length(T.idx)
  
  S1=as.numeric(t(diff.i)%*%Psi.ii.inv%*%diff.i)
  fact.cov=(nu.cur+S1-2)/df.o ## (nu+S1-2)/(df.o-2)*(df.o-2)/df.o
  ## factor for scale matrix * factor for prediction
  
  pred.scale.mtx=fact.cov*tilde.Psi.oo
  pred.scale.mtx=(pred.scale.mtx+t(pred.scale.mtx))/2
  
  Z.sample.out.new=rmvt(1,sigma=pred.scale.mtx,df=df.o,delta=tilde.mu.o,
                        type='shifted')
  Z.sample.new[-idx.sub]=Z.sample.out.new
  Z.sample.new[idx.sub]=Z.sample.in.new
  
  prob.sample.new=expit(Z.sample.new+rnorm(length(Z.sample.new),mean=0,sd=sqrt(sigma.epsilon.sample)))
  
  return(prob.sample.new)
}
