#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <limits>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
vec plog(vec p){
  for (unsigned i = 0;i < p.n_elem; i++){
    if(is_finite(log(p[i]))&is_finite(log(1-p[i]))){
      p[i]=p[i]*log(p[i])+(1-p[i])*log(1-p[i]);
    } else p[i] = 0;
  }
  return p;
}

// [[Rcpp::export]]
List Estimate(const vec & residual,const mat & X,vec m,vec s2,vec p,double s2_eps,double s2_b, double pi){
  const unsigned int n = X.n_rows;
  const unsigned int r = X.n_cols;
  mat L = zeros<mat>(101,100); 
  
  const double res2 = dot(residual,residual);
  const vec XTX = sum(X%X).t();
  vec mp = m%p;
  vec mpmp = mp%mp;
  vec Xmp = X*mp;
  vec m2s2p = (m%m+s2)%p;
  
  for(unsigned int iter_outer = 0;iter_outer < L.n_cols;iter_outer++){
    double Ey_Xb2 = res2-2*dot(residual,Xmp)+dot(Xmp,Xmp)+dot(XTX,m2s2p-mpmp);
    double Eq_likelihood = -(n*log(s2_eps)+Ey_Xb2/s2_eps)/2;
    double Eq_prior = -sum(m2s2p)/(2*s2_b)+sum(p)*(-log(s2_b)/2-log(2*datum::pi)/2+log(pi)-log(1-pi))+r*log(1-pi);
    double Eq_q = sum(-p%log(s2)/2-p/2+plog(p));

    s2 = 1/(XTX/s2_eps+1/s2_b);
    m2s2p = (m%m+s2)%p;
    L(0,iter_outer) = Eq_likelihood+Eq_prior-Eq_q;
    for(unsigned int iter_inner = 1;iter_inner < L.n_rows-1;iter_inner++){
      for (unsigned int k = 0; k < r;k++){
        Xmp -= mp[k]*X.col(k); 
        m[k] = s2[k]*(dot(X.col(k),residual-Xmp))/s2_eps; 
        p[k] = 1-1/(1+sqrt(s2[k]/s2_b)*pi/(1-pi)*exp(m[k]*m[k]/(2*s2[k])));
        mp[k] = m[k] * p[k];
        Xmp += mp[k]*X.col(k);
      } 
      
      mpmp = mp%mp;
      m2s2p = (m%m+s2)%p;
      Ey_Xb2 = res2-2*dot(residual,Xmp)+dot(Xmp,Xmp)+dot(XTX,m2s2p-mpmp);
      Eq_likelihood = -(n*log(s2_eps)+Ey_Xb2/s2_eps)/2;
      Eq_prior = -sum(m2s2p)/(2*s2_b)+sum(p)*(-log(s2_b)/2-log(2*datum::pi)/2+log(pi)-log(1-pi))+r*log(1-pi);
      Eq_q = sum(-p%log(s2)/2-p/2+plog(p));
      L(iter_inner,iter_outer) = Eq_likelihood+Eq_prior-Eq_q;
      if((L(iter_inner,iter_outer)-L(iter_inner-1,iter_outer))<0.000001) break;
    } 
    s2_eps=Ey_Xb2/n;
    pi=mean(p);
    s2_b=sum(m2s2p)/(r*pi);
    
    Eq_likelihood = -(n*log(s2_eps)+Ey_Xb2/s2_eps)/2;
    Eq_prior = -sum(m2s2p)/(2*s2_b)+sum(p)*(-log(s2_b)/2-log(2*datum::pi)/2+log(pi)-log(1-pi))+r*log(1-pi);
    L(L.n_rows-1,iter_outer) = Eq_likelihood+Eq_prior;
    if(iter_outer > 0) if((L(L.n_rows-1,iter_outer)-L(L.n_rows-1,iter_outer-1))<0.000001) break;
  }
  
  return List::create(_("L") = L, _("m") = m , _("s2") = s2, _("p") = p, _("s2_eps") = s2_eps, _("s2_b") = s2_b, _("pi") = pi);
}

/*** R
gc()
while(!require('rstudioapi')) install.packages('rstudioapi')
while(!require('xtable')) install.packages('xtable')
while(!require('ggplot2')) install.packages('ggplot2')
while(!require('reshape2')) install.packages('reshape2')
while(!require('mvtnorm')) install.packages('mvtnorm')
while(!require('dplyr')) install.packages('dplyr')
while(!require('glmnet')) install.packages('glmnet')
setwd(dirname(getSourceEditorContext()$path))

load('GWAS.Rdata')
n=nrow(L$X)
p=ncol(L$X)

Result=vector('list',4)
for(i in 1:ncol(L$y)){
  start_time=Sys.time()
  parameter=Estimate(L$y[,i],L$X,rep(0,p),rep(1,p),rep(0.5,p),5,1,min(n*0.1/p,0.5))
  end_time=Sys.time()
  parameter$outer_L=parameter$L[nrow(parameter$L),parameter$L[nrow(parameter$L),]!=0]
  parameter$L=apply(parameter$L[-nrow(parameter$L),],2,function(x) x[x!=0])
  parameter$L=parameter$L[sapply(parameter$L,length)>0]
  plot(parameter$outer_L,type='l',ylab='Likelihood',xlab='Step')
  plot(log(parameter$p),ylab=expression('log'[10]*'(posterior probability)'))
  print(end_time-start_time)
  Result[[i]]=list(parameter,start_time,end_time,index = which(parameter$p>0.5))
  save(parameter,start_time,end_time,index,n,p,file=paste0('GWAS',i,'.Rdata'))
}

# Table
lapply(Result,function(x) data.frame(x$index,x[[1]]$m[x$index],x[[1]]$s2[x$index],x[[1]]$p[x$index]))%>%bind_rows%>%xtable(digits=6)
lapply(Result,function(x) data.frame('b'=x[[1]]$s2_b,'e'=x[[1]]$s2_eps,'p'=p*x[[1]]$pi))%>%bind_rows%>%xtable(digits = 6)


# Plot 
lapply(Result,function(x) 1:length(x[[1]]$outer_L))%>%melt%>%cbind(likelihood=lapply(Result,function(x) x[[1]]$outer_L)%>%unlist)%>%
  ggplot(aes(x=value,y=likelihood,color=as.factor(L1)))+
  geom_line()+
  labs(color="Phenotype")+
  xlab("Step")+
  ylab("Likelihood")
sapply(Result,function(x) log(x[[1]]$p))%>%melt%>%
  ggplot(aes(x=Var1,y=value))+
  geom_point()+
  facet_grid(Var2~.)+
  xlab("Index")+
  ylab(expression('log'[10]*'(posterior probability)'))
LASSO=vector('list',4)

for(i in 1:ncol(L$y)){
  gc()
  LASSO[[i]]=cv.glmnet(L$X,L$y[,i],intercept = FALSE)
}

# Result LASSO
sapply(LASSO,function(z) c(which.min(z$cvm)%>%
                             (function(x) c(z$lambda[x],z$cvm[x],z$cvsd[x],z$nzero[x])),
                           which.max(z$lambda*(z$cvm<z$cvup[which.min(z$cvm)]))%>%
                             (function(x) c(z$lambda[x],z$cvm[x],z$cvsd[x],z$nzero[x])))
)%>%t%>%xtable(digits = 6)

# Table
data.frame(LASSO=sapply(LASSO,function(z) z$cvm[which.max(z$lambda*(z$cvm<z$cvup[which.min(z$cvm)]))]),
           VEM=colMeans((y-sapply(Result,function(x) X[,x$index]%*%x[[1]]$m[x$index]))^2))%>%xtable(digits=6)
*/

