## Script for performing modeling the effect random then greedy boosting

library(tidyverse)
library(glmnet)
library(mvtnorm)


## Implementing random then greedy approach to FS-epsilon

 random_then_greedyFSe<-function(X,y,Xte,yte,Iters,subset.size=20,step.size=0.01) {
  
  n=nrow(X); p=ncol(X);
  
  subset.size<-min(p,subset.size);
  
  betak<-rep(0,p);
  sparsity<-train.errs<-test.errs<-rep(0,ceiling(Iters/50));
  
  counter<-0;
  res<-y;
  
  for (i in 1:Iters){
    
    ids<-sample(seq(1,p),subset.size,replace=FALSE)  

    grad<- -t(X[,ids])%*%res; ik<-which.max(abs(grad));
    
    if (i/50==ceiling(i/50)) {  
      counter<-counter+1;
      train.errs[counter]<- sum((y-X%*%betak)^2); test.errs[counter]<-sum((yte-Xte%*%betak)^2);
      sparsity[counter]<-sum(abs(betak) > 1e-6); 
                             }
    
    betak[ids[ik]] <- betak[ids[ik]] - step.size*sign(grad[ik])
    res<- res + step.size*sign(grad[ik])*X[,ids[ik]];
    
                  }
  
  train.errs<-train.errs[1:counter];test.errs<-test.errs[1:counter]; sparsity<-sparsity[1:counter];
  
  return(list(train.errs=train.errs,test.errs=test.errs,betak=betak,sparsity=sparsity))
  
 }
 
 generate_cov <- function(type="ar1", p, rho){
   ## Generate covariate matrix
   if(type == "ar1"){
     Sigma = matrix(0,nc=p,nr=p)
     for(i in 1:p){
       for(j in 1:p){
         Sigma[i,j] = rho^(abs(i-j))
       }
     }
   } else if (type == "exchengable"){
     Sigma = matrix(rho, nc=p,nr=p)
     diag(Sigma) = 1
   } else if(type == "identity"){
     Sigma = diag(p)
   }
   return(Sigma)
 }
 


## To perform simulations with different values of rho
 n = 500; p = 5000
###################################################
## Exchengable covariance matrix
###################################################
#options(repr.plot.width=10, repr.plot.height=10)
rho_vals = c(0,0.5,0.99)
par(mfrow=c(length(rho_vals),2))
for (r in 1:length(rho_vals)){
    rho = rho_vals[r]
    message("rho:",rho,"\n")
    ## Generate data
    Sigma = generate_cov(type="exchengable",p, rho=rho)
    s.percent = 0.20
    X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)
    beta = 2*rnorm(p)
    loc.signal <-sample.int(p,floor(s.percent*p))
    beta[setdiff(1:p, loc.signal)] = 0
    Y = X %*% beta
    
    ## Partition
    ## partition into testing and training sets
    tr_id = seq(1, floor(0.80*n))
    test_id = setdiff(1:n, tr_id)

    X_train = X[tr_id,]; Y_train = Y[tr_id]
    X_test = X[test_id,]; Y_test = Y[test_id] 

    ## perform scaling
    m1<-apply(X_train,2,mean); s1<-apply(X_train,2,sd)*sqrt(length(tr_id));
    X_train<-scale(X_train,center=m1,scale=s1)
    X_test <-scale(X_test,center=m1,scale=s1)
    
    ## Run random then greedy FSe

    ## RtG
    Iters=50000
    ptm <- proc.time()
    out1<-random_then_greedyFSe(X_train,Y_train,X_test,Y_test,Iters,subset.size=20,step.size=0.1)
    ptm1<-proc.time()-ptm
    time1<-seq(1,length(out1$train.errs))*ptm1[1]/length(out1$train.errs);

    ## Random
    ptm <- proc.time()
    out2<-random_then_greedyFSe(X_train,Y_train,X_test,Y_test,Iters,subset.size=1,step.size=0.1)
    ptm2<-proc.time()-ptm
    time2<-seq(1,length(out2$train.errs))*ptm2[1]/length(out2$train.errs);

    ## Greedy
    ptm <- proc.time()
    out3<-random_then_greedyFSe(X_train,Y_train,X_test,Y_test,Iters,subset.size=p,step.size=0.1)
    ptm3<-proc.time()-ptm
    time3<-seq(1,length(out3$train.errs))*ptm3[1]/length(out3$train.errs);
    
    ## Plotting per iters
    plot(out3$train.errs, col="red",xlab="Iterations",ylab="Training Error",main=paste0("Train vs Iters, rho=",rho));
    points(out2$train.errs, col="blue");
    points(out1$train.errs, col="orange")
    legend("topright", legend=c("subset=p","subset=1","subset=20"),col=c("red","blue","orange"),lty=1,pch=19)

    plot(out3$test.errs, col="red",xlab="Iterations",ylab="Test Error",main=paste0("Test vs Iters, rho=",rho));
    points(out2$test.errs, col="blue");
    points(out1$test.errs, col="orange")
    legend("topright", legend=c("subset=p","subset=1","subset=20"),col=c("red","blue","orange"),lty=1,pch=19)

}


###############################
## AR(1) covariance 
###############################

## To perform simulations with different values of rho
#options(repr.plot.width=10, repr.plot.height=10)
rho_vals = c(0,0.5,0.99)
par(mfrow=c(length(rho_vals),2))
for (r in 1:length(rho_vals)){
    rho = rho_vals[r]
    message("rho:",rho,"\n")
    ## Generate data
    Sigma = generate_cov(type="ar1",p, rho=rho)
    s.percent = 0.20
    X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)
    beta = 2*rnorm(p)
    loc.signal <-sample.int(p,floor(s.percent*p))
    beta[setdiff(1:p, loc.signal)] = 0
    Y = X %*% beta
    
    ## Partition
    ## partition into testing and training sets
    tr_id = seq(1, floor(0.80*n))
    test_id = setdiff(1:n, tr_id)

    X_train = X[tr_id,]; Y_train = Y[tr_id]
    X_test = X[test_id,]; Y_test = Y[test_id] 

    ## perform scaling
    m1<-apply(X_train,2,mean); s1<-apply(X_train,2,sd)*sqrt(length(tr_id));
    X_train<-scale(X_train,center=m1,scale=s1)
    X_test <-scale(X_test,center=m1,scale=s1)
    
    ## Run random then greedy FSe

    ## RtG
    Iters=50000
    ptm <- proc.time()
    out1<-random_then_greedyFSe(X_train,Y_train,X_test,Y_test,Iters,subset.size=20,step.size=0.1)
    ptm1<-proc.time()-ptm
    time1<-seq(1,length(out1$train.errs))*ptm1[1]/length(out1$train.errs);

    ## Random
    ptm <- proc.time()
    out2<-random_then_greedyFSe(X_train,Y_train,X_test,Y_test,Iters,subset.size=1,step.size=0.1)
    ptm2<-proc.time()-ptm
    time2<-seq(1,length(out2$train.errs))*ptm2[1]/length(out2$train.errs);

    ## Greedy
    ptm <- proc.time()
    out3<-random_then_greedyFSe(X_train,Y_train,X_test,Y_test,Iters,subset.size=p,step.size=0.1)
    ptm3<-proc.time()-ptm
    time3<-seq(1,length(out3$train.errs))*ptm3[1]/length(out3$train.errs);
    
    ## Plotting per iters
    plot(out3$train.errs, col="red",xlab="Iterations",ylab="Training Error",main=paste0("Train vs Iters, rho=",rho));
    points(out2$train.errs, col="blue");
    points(out1$train.errs, col="orange")
    legend("topright", legend=c("subset=p","subset=1","subset=20"),col=c("red","blue","orange"),lty=1,pch=19)

    plot(out3$test.errs, col="red",xlab="Iterations",ylab="Test Error",main=paste0("Test vs Iters, rho=",rho));
    points(out2$test.errs, col="blue");
    points(out1$test.errs, col="orange")
    legend("topright", legend=c("subset=p","subset=1","subset=20"),col=c("red","blue","orange"),lty=1,pch=19)

}



