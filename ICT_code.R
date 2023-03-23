
#########################################################################
#
# This file contains R functions implementing the estimators described  
# in the paper: Kowalczyk B., Niemiro W., Wieczorkowski R.:

#   "Item count technique based on two different treatment groups"
#   (Item count technique with continuous auxiliary variable)
# (submitted to Statistical Methods in Medical Research)
#
# actualization date: 02.05.2017
#


### section for X - normal distribution

estep.norm <- function(y1,y2,pp,mi,sigma2,a=1)
#
# Expectation step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   pp - estimate for parameter pi (sensitive fraction)
#   mi - estimate for parameter mi
#   sigma2 - estimate for parameter sigma2
#   a - additional parameter of the method (multiplier for Z)
#
# Values:
#   output list contains two vectors with expected values 
#
{
  z1_estep <- pp/(pp+(1-pp)*exp(a*(2*y1-2*mi+a)/(2*sigma2)))
  z2_estep <- pp/(pp+(1-pp)*exp(a*(-2*y2+2*mi+a)/(2*sigma2)))
  #return(cbind(z1_estep,z2_estep))
  return(list(z1_estep,z2_estep) )
}


mstep.norm <- function(y1,y2,e.step,a=1)
#
# Maximization step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   e.step - values from Expectation step
#   a - additional parameter of the method (multiplier for Z)
#
# Values:
#   output list contains values of parameters (pi,mi,sigma2) 
#   from M-step of EM algorithm
{
  n1<-length(y1)
  n2<-length(y2)
  
  #z1<-e.step[,1]
  #z2<-e.step[,2]
  z1<-e.step[[1]]
  z2<-e.step[[2]]
  
  # estimate pi
  pi_temp <- (sum(z1)+sum(z2))/(n1+n2)
  
  # estimate mi
  mi_temp <-(sum(y1+a*z1)+sum(y2-a*z2))/(n1+n2) 
  
  # estimate sigma2
  #sigma2_temp <- (sum((y1+z1-mi_temp)^2) + sum((y1-z2-mi_temp)^2))/(n1+n2) #old
  sigma2_temp <- sum(z1*((y1+a-mi_temp)^2) + (1-z1)*((y1-mi_temp)^2)) +
                 sum(z2*((y2-a-mi_temp)^2) + (1-z2)*((y2-mi_temp)^2))
  sigma2_temp <- sigma2_temp/(n1+n2)
  #cat("sigma2_temp = ",sigma2_temp,"\n")
  
  list(pi_temp,mi_temp,sigma2_temp)   
}



em.algo.norm <- function(y1,y2,pi_inits=NULL,mi_inits=NULL,sigma2_inits=NULL,maxit=10000,
                         tol=1e-6,a=1,info=FALSE,AIC=FALSE)
#
# EM algorithm for Maximum likelihood estimation
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   pi_inits - initial value for probability pi
#   mi_inits - initial vlues for parameter mi
#   sigma2_inits - initial vlues for parameter sigma2
#   a - additional parameter of the method (multiplier for Z)
#   maxit - maximal number of iterations (default 10000)
#   tol - values for testing convergence of algorithm (default 1-6)
#   info - logical value for printing information about iterations
#   AIC - if TRUE then AIC measure will be computed
#
# Values:
#   output list contains ML estimators of parameters (pi,mi,sigma2) 
# 
{
  # Initial parameter estimates
  flag <- 0
  
  if (is.null(pi_inits))
  {
    if (info) cat("Initialization of parameters by method of moments ...\n")
    n1<-length(y1)
    n2<-length(y2)
    pi_inits<-(mean(y2)-mean(y1))/(2*a)
    pi_inits <- ifelse(pi_inits<0 | pi_inits>1,0.5,pi_inits)
    mi_inits<-(mean(y2)+mean(y1))/2
    sigma2_inits<-(var(y1)*(n1-1)/n1+var(y2)*(n2-1)/n2)/2 - a*a*pi_inits*(1-pi_inits)
    if (info) cat("pi_inits=",pi_inits," mi_inits=",mi_inits,"sigma2_inits=",sigma2_inits,"\n")
  }
  pi_cur <- pi_inits; mi_cur <- mi_inits; sigma2_cur <- sigma2_inits;
  
  
  # Iterate between expectation and maximization steps
  for(i in 1:maxit){
    cur <- c(pi_cur,mi_cur,sigma2_cur)
    #print(cur)
    new <- mstep.norm(y1,y2,estep.norm(y1,y2, pi_cur, mi_cur, sigma2_cur,a=a),a=a)
    pi_new <- new[[1]]; mi_new <- new[[2]]; sigma2_new <- new[[3]];
	if (is.na(pi_new)) pi_new<-pi_cur
	if (is.na(mi_new)) mi_new<-mi_cur
	if (is.na(sigma2_new)) sigma2_new<-sigma2_cur
	
    new_step <- c(pi_new,mi_new,sigma2_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    #cat("blad = ",round(abs(cur - new_step),4),"\n")    
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}
    
    
    # Otherwise continue iteration
    pi_cur <- pi_new; mi_cur <- mi_new; sigma2_cur <- sigma2_new;
  }
  if(!flag) warning("Didn’t converge\n")
  if (info) cat(i," iterations \n")
  
  if (AIC) 
  {
    fL<-function(x)
    {
      return(log(dnorm(x,mean=mi_cur,sd=sqrt(sigma2_cur))))
    }
  
    z<-estep.norm(y1,y2,pi_cur,mi_cur,sigma2_cur,a=a)
    z1<-z[[1]]
    z2<-z[[2]]
  
    lgL<-log(pi_cur)*(sum(z1)+sum(z2))+log(1-pi_cur)*(sum(1-z1)+sum(1-z2))+
      sum(z1*fL(y1+a))+sum((1-z1)*fL(y1))+sum(z2*fL(y2-a))+sum((1-z2)*fL(y2))

    lgL<-lgL-sum(z1*log(z1))-sum((1-z1)*log(1-z1))-sum(z2*log(z2))-sum((1-z2)*log(1-z2))
      
    k_param<-2
    AICvalue<-2*(k_param+1)-2*lgL
    out<-list(pi=pi_cur,mi=mi_cur,sigma2=sigma2_cur,AIC=AICvalue)
  }
  else out<-list(pi=pi_cur, mi=mi_cur, sigma2=sigma2_cur)
  
  return(out)
}









### section for X - lognormal distribution

estep.lognorm <- function(y1,y2,pp,mi,sigma2,a=1)
#
# Expectation step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   pp - estimate for parameter pi (sensitive fraction)
#   mi - estimate for parameter mi
#   sigma2 - estimate for parameter sigma2
#   a - additional parameter of the method (multiplier for Z)
#
# Values:
#    output list contains two vectors with expected values 
#
{
  z1_estep <- (pp/(y1+a))/
    ( 
      pp/(y1+a)+((1-pp)/(y1))*exp((log(y1+a)-log(y1))*(log(y1+a)+log(y1)-2*mi)/(2*sigma2)) 
    )
  z2_estep <- (pp/(y2-a))/
    ( 
      pp/(y2-a)+((1-pp)/(y2))*exp((log(y2-a)-log(y2))*(log(y2-a)+log(y2)-2*mi)/(2*sigma2)) 
    )
  
  z1_estep <- ifelse(y1<=0,1,z1_estep)
  z2_estep <- ifelse(y2<=a,0,z2_estep)
  
  return(list(z1_estep,z2_estep))
}



mstep.lognorm <- function(y1,y2,e.step,a=1)
#
# Maximization step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   e.step - values from Expectation step
#   a - additional parameter of the method (multiplier for Z)
#
# Values:
#   output list contains values of parameters (pi,mi,sigma2) 
#   from M-step of EM algorithm
{
  n1<-length(y1)
  n2<-length(y2)
  
  z1<-e.step[[1]]
  z2<-e.step[[2]]
  
  # estimate pi
  pi_temp <- (sum(z1)+sum(z2))/(n1+n2)
  
  # estimate mi
  mi_temp <-( sum(z1*log(y1+a)+(1-z1)*log(y1)) + sum(z2*log(y2-a)+(1-z2)*log(y2)) )/(n1+n2) 
  
  # estimate sigma2
  #sigma2_temp <- (sum((y1+z1-mi_temp)^2) + sum((y1-z2-mi_temp)^2))/(n1+n2) #old
  sigma2_temp <- sum(z1*((log(y1+a)-mi_temp)^2) + (1-z1)*((log(y1)-mi_temp)^2)) +
    sum(z2*((log(y2-a)-mi_temp)^2) + (1-z2)*((log(y2)-mi_temp)^2))
  sigma2_temp <- sigma2_temp/(n1+n2)
  #cat("sigma2_temp = ",sigma2_temp,"\n")
  
  list(pi_temp,mi_temp,sigma2_temp)   
}



em.algo.lognorm <- function(y1,y2,pi_inits=NULL,mi_inits=NULL,sigma2_inits=NULL,maxit=10000,tol=1e-6,a=1,info=FALSE,AIC=FALSE)
#
# EM algorithm for Maximum likelihood estimation
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   pi_inits - initial value for probability pi
#   mi_inits - initial vlues for parameter mi
#   sigma2_inits - initial vlues for parameter sigma2
#   a - additional parameter of the method (multiplier for Z)
#   maxit - maximal number of iterations (default 10000)
#   tol - values for testing convergence of algorithm (default 1-6)
#   info - logical value for printing information about iterations
#   AIC - if TRUE then AIC measure will be computed
#
# Values:
#   output list contains ML estimators of parameters (pi,mi,sigma2) 
# 
{
  # Initial parameter estimates
  flag <- 0
  if (is.null(pi_inits))
  {
    if (info) cat("Initialization of parameters by method of moments ...\n")
    n1<-length(y1)
    n2<-length(y2)
    pi_inits<-(mean(y2)-mean(y1))/(2*a)
    pi_inits <- ifelse(pi_inits<0 | pi_inits>1,0.5,pi_inits)
  
    B<-(mean(y2)+mean(y1))/2
    A<-(var(y1)*(n1-1)/n1+var(y2)*(n2-1)/n2)/2 - a*a*pi_inits*(1-pi_inits)
    mi_inits<-log(B)-0.5*log(A/B^2+1)
    sigma2_inits<-log(A/B^2+1)
    if (info) cat("pi_inits=",pi_inits," mi_inits=",mi_inits,"sigma2_inits=",sigma2_inits,"\n")
  }
  
  pi_cur <- pi_inits; mi_cur <- mi_inits; sigma2_cur <- sigma2_inits;
  
  # Iterate between expectation and maximization steps
  for(i in 1:maxit){
    cur <- c(pi_cur,mi_cur,sigma2_cur)
    #print(cur)
    new <- mstep.lognorm(y1,y2,estep.lognorm(y1,y2, pi_cur, mi_cur, sigma2_cur,a=a),a=a)
    pi_new <- new[[1]]; mi_new <- new[[2]]; sigma2_new <- new[[3]];
    new_step <- c(pi_new,mi_new,sigma2_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    #cat("blad = ",round(abs(cur - new_step),4),"\n")    
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}
    
    
    # Otherwise continue iteration
    pi_cur <- pi_new; mi_cur <- mi_new; sigma2_cur <- sigma2_new;
  }
  if(!flag) warning("Didn’t converge\n")
  if (info) cat(i," iterations \n")

  if (AIC) 
  {
    fL<-function(x)
    {
      w<-log(dlnorm(x,meanlog=mi_cur,sdlog=sqrt(sigma2_cur)))
      w[is.nan(w)]<-0
      return(w)
    }
    
    z<-estep.lognorm(y1,y2,pi_cur,mi_cur,sigma2_cur,a=a)
    z1<-z[[1]]
    z2<-z[[2]]
    
    lgL<-log(pi_cur)*(sum(z1)+sum(z2))+log(1-pi_cur)*(sum(1-z1)+sum(1-z2))+
      sum(z1*fL(y1+a))+sum((1-z1)*fL(y1))+sum(z2*fL(y2-a))+sum((1-z2)*fL(y2))
    
    lgL<-lgL-sum(z1*log(z1))-sum((1-z1)*log(1-z1))-sum(z2*log(z2))-sum((1-z2)*log(1-z2))
    
    k_param<-2
    AICvalue<-2*(k_param+1)-2*lgL
    out<-list(pi=pi_cur,mi=mi_cur,sigma2=sigma2_cur,AIC=AICvalue)
  }
  else out=list(pi=pi_cur, mi=mi_cur, sigma2=sigma2_cur)
  
  return(out)
}




### section for X - gamma distribution

estep.gamma <- function(y1,y2,pp,k,theta,a=1)
#
# Expectation step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   pp - estimate for parameter pi (sensitive fraction)
#   k - estimate for parameter k
#   theta - estimate for parameter theta
#   a - additional parameter of the method (multiplier for Z)
#
# Values:
#    output list contains two vectors with expected values 
#
{
  z1_estep <- (pp*((y1+a)^(k-1)))/( pp*((y1+a)^(k-1)) + (1-pp)*exp(a/theta)*(y1^(k-1) ) )
  z1_estep <- ifelse(y1<=0,1,z1_estep)
  z2_estep <- (pp*exp(a/theta)*((y2-a)^(k-1)))/( pp*exp(a/theta)*((y2-a)^(k-1)) + (1-pp)*(y2^(k-1) ) )
  z2_estep <- ifelse(y2<=a,0,z2_estep)
  return(list(z1_estep,z2_estep))
}


mstep.gamma <- function(y1,y2,e.step,B,S1,S2,a=1)
#
# Maximization step for Maximum likelihood estimation using EM algorithm
#
# Arguments:
#   y1 - sample from first treatment group (X-Z)
#   y2 - sample from second treatment group (X+Z)
#   e.step - values from Expectation step
#   B, S1, S2 - values computed only once in initial step of EM 
#   a - parameter
#
# Values:
#   output list contains values of parameters (pi,k,theta) 
#   from M-step of EM algorithm
{
  n1<-length(y1)
  n2<-length(y2)
  
  z1<-e.step[[1]]
  z2<-e.step[[2]]
  
  # estimate pi
  pi_temp <- (sum(z1)+sum(z2))/(n1+n2)
  
  
  # estimate k (B, S1, S2 - depend on y !)
  #B <- (mean(y1)+mean(y2))/2
  #S1 <- sum((y1-mean(y1))^2)/n1
  #S2 <- sum((y2-mean(y2))^2)/n2
  A <- 0.5*(S1+S2) - a*a*pi_temp*(1-pi_temp)
  k_temp <- B*B/A
  
  # estimate theta
  theta_temp <-(sum(y1+a*z1)+sum(y2-a*z2))/(k_temp*(n1+n2)) 
  
  list(pi_temp,k_temp,theta_temp)   
}


em.algo.gamma <- function(y1,y2,pi_inits=NULL,k_inits=NULL,theta_inits=NULL,a=1,maxit=10000,tol=1e-6,info=FALSE,AIC=FALSE)
#
# EM algorithm for Maximum likelihood estimation
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   pi_inits - initial value for probability pi
#   k_inits - initial vlues for parameter k
#   theta_inits - initial vlues for parameter theta
#   a - additional parameter of the method (multiplier for Z)
#   maxit - maximal number of iterations (default 10000)
#   tol - values for testing convergence of algorithm (default 1-6)
#   info - logical value for printing information about iterations
#   AIC - if TRUE then AIC measure will be computed
#
# Values:
#   output list contains ML estimators of parameters (pi,k,theta) 
# 
{
  # Initial parameter estimates
  flag <- 0
  if (is.null(pi_inits))
  {
    if (info) cat("Initialization of parameters by method of moments ...\n")
    n1<-length(y1)
    n2<-length(y2)
    pi_inits<-(mean(y2)-mean(y1))/(2*a)
    pi_inits <- ifelse(pi_inits<0 | pi_inits>1,0.5,pi_inits)
    
    B <- (mean(y1)+mean(y2))/2
    S1 <- sum((y1-mean(y1))^2)/n1
    S2 <- sum((y2-mean(y2))^2)/n2
    A <- 0.5*(S1+S2) - a*a*pi_inits*(1-pi_inits)
    k_inits <- B*B/A
    theta_inits<-A/B
    if (info) cat("pi_inits=",pi_inits," k_inits=",k_inits,"theta_inits=",theta_inits,"\n")
  }
  
  pi_cur <- pi_inits; k_cur <- k_inits; theta_cur <- theta_inits;
  
  B <- (mean(y1)+mean(y2))/2
  S1 <- sum((y1-mean(y1))^2)/length(y1)
  S2 <- sum((y2-mean(y2))^2)/length(y2)
  
  # Iterate between expectation and maximization steps
  for(i in 1:maxit){
    cur <- c(pi_cur,k_cur,theta_cur)
    new <- mstep.gamma(y1,y2,estep.gamma(y1,y2, pi_cur, k_cur, theta_cur,a=a),B,S1,S2,a=a)
    pi_new <- new[[1]]; k_new <- new[[2]]; theta_new <- new[[3]];
    new_step <- c(pi_new,k_new,theta_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}
     
    # Otherwise continue iteration
    pi_cur <- pi_new; k_cur <- k_new; theta_cur <- theta_new;
  }
  if(!flag) warning("Didn’t converge\n")
  if (info) cat(i," iterations \n")

  if (AIC) 
  {
    fL<-function(x)
    {
      w<-log(dgamma(x,shape=k_cur,scale=theta_cur))
      w[is.nan(w)]<-0
      return(w)
    }
    
    z<-estep.gamma(y1,y2,pi_cur,k_cur,theta_cur,a=a)
    z1<-z[[1]]
    z2<-z[[2]]
    
    lgL<-log(pi_cur)*(sum(z1)+sum(z2))+log(1-pi_cur)*(sum(1-z1)+sum(1-z2))+
      sum(z1*fL(y1+a))+sum((1-z1)*fL(y1))+sum(z2*fL(y2-a))+sum((1-z2)*fL(y2))
    
    lgL<-lgL-sum(z1*log(z1))-sum((1-z1)*log(1-z1))-sum(z2*log(z2))-sum((1-z2)*log(1-z2))
    
    k_param<-2
    AICvalue<-2*(k_param+1)-2*lgL
    out<-list(pi=pi_cur,k=k_cur,theta=theta_cur,AIC=AICvalue)
  }
  else out<-list(pi=pi_cur, k=k_cur, theta=theta_cur)
  
  return(out)
}




#######

### section for X - Poisson distribution

estep.poi <- function(y1,y2,pp,lambda,a=1)
  #
  # Expectation step for Maximum likelihood estimation using EM algorithm
  #
  # Arguments:
  #   y1 - sample from first treatment group (X-Z)
  #   y2 - sample from second treatment group (X+Z)
  #   pp - estimate for parameter pi
  #   lambda - estimate for parameter lambda
  #   a - parameter
  #
  # Values:
#   output data frame contains two columns with expected values 
#

{
  z1_estep <- (pp*lambda^a)/(pp*lambda^a+(factorial(y1+a)/factorial(y1))*(1-pp))
  z1_estep <- ifelse(y1<0,1,z1_estep)
  z2_estep <- (pp*factorial(y2)/factorial(y2-a))/(pp*factorial(y2)/factorial(y2-a)+(lambda^a)*(1-pp))
  z2_estep <- ifelse(y2<a,0,z2_estep)
  
  return(list(z1_estep,z2_estep))
}


mstep.poi <- function(y1,y2,e.step,a=1)
  #
  # Maximization step for Maximum likelihood estimation using EM algorithm
  #
  # Arguments:
  #   y1 - sample from first treatment group (X-Z)
  #   y2 - sample from second treatment group (X+Z)
  #   e.step - values from Expectation step
  #   a - parameter
  #
  # Values:
  #   output list contains values of parameters (pi,lambda) 
#   from M-step of EM algorithm
{
  
  # estimate pi
  pi_temp <- mean(unlist(e.step))
  
  # estimate lambda
  lambda_temp <-(sum(y1+a*e.step[[1]])+sum(y2-a*e.step[[2]]))/(length(y1)+length(y2)) 
  
  list(pi_temp,lambda_temp)   
}



em.algo.poi <- function(y1,y2,pi_inits=NULL,lambda_inits=NULL,a=1,maxit=10000,tol=1e-6,info=FALSE)
#
# EM algorithm for Maximum likelihood estimation
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   pi_inits - initial value for probability pi
#   lambda_inits - initial vlues for parameter lambda
#   a - additional parameter of the method (multiplier for Z)
#   maxit - maximal number of iterations (default 10000)
#   tol - values for testing convergence of algorithm (default 1-6)
#   info - logical value for printing information about iterations
#
# Values:
#   output list contains ML estimators of parameters (pi,lambda) 
# 
{
  # Initial parameter estimates
  flag <- 0
  
  if (is.null(pi_inits))
  {
    if (info) cat("Initialization of parameters by method of moments ...\n")
    pi_inits<-(mean(y2)-mean(y1))/(2*a)
    pi_inits <- ifelse(pi_inits<0 | pi_inits>1,0.5,pi_inits)
    lambda_inits<-(mean(y2)+mean(y1))/2
    if (info) cat("pi_inits=",pi_inits," lambda_inits=",lambda_inits,"\n")
  }
  
  pi_cur <- pi_inits; lambda_cur <- lambda_inits;
  
  # Iterate between expectation and maximization steps
  for(i in 1:maxit){
    cur <- c(pi_cur,lambda_cur)
    new <- mstep.poi(y1,y2,estep.poi(y1,y2, pi_cur, lambda_cur,a=a),a=a)
    #new <- mstep.poi1(y1,y2,estep.poi(y1,y2, pi_cur, lambda_cur))
    pi_new <- new[[1]]; lambda_new <- new[[2]];
    
    #if (is.na(pi_new)) pi_new<-pi_cur
    #if (is.na(lambda_new)) lambda_new<-lambda_cur
    
    new_step <- c(pi_new,lambda_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}
    
    
    # Otherwise continue iteration
    pi_cur <- pi_new; lambda_cur <- lambda_new;
  }
  if(!flag) warning("Didn’t converge\n")
  if (info) cat(i," iterations \n")
  
  return(list(pi=pi_cur, lambda=lambda_cur))
}






ML.boot.ci<-function(y1,y2,a=1,conf=0.95,model="normal")
#
# Bootstrap nonparametric confidence intervals for Maximum likelihood estimation
# of proportion 'pi', for 3 models
#
# Arguments:
#   y1 - sample from first treatment group (X-a*Z)
#   y2 - sample from second treatment group (X+a*Z)
#   a - additional parameter of the method (multiplier for Z)
#   conf - level of confidence intervals
#   model - type of distribution for X: normal, lognormal or gamma
#
# Values:
#   output list contains lower and upper values for confidence interval 
# 
{
  require(boot)
  n1<-length(y1)
  n2<-length(y2)
  
  boot.stat<-function(data,indices,a=a)
  {
    #xx<-data[indices,1]
    #yy<-data[indices,2]
    xx<-data[indices[1:n1]]
    yy<-data[indices[(n1+1):(n1+n2)]]
    
    if (model=="normal") wy<-em.algo.norm(xx,yy,a=a)
    else if (model=="lognormal") wy<-em.algo.lognorm(xx,yy,a=a)
    else if (model=="gamma") wy<-em.algo.gamma(xx,yy,a=a)
    else if (model=="Poisson") wy<-em.algo.poi(xx,yy,a=a)
    else { cat("Bad name of the model distribution ...\n"); stop }
    
    return(as.vector(wy[[1]]))
  }
  
  #d<-cbind(y1,y2)
  d<-c(y1,y2)
  b<-boot(data=d,statistic=boot.stat,strata=c(rep(1,n1),rep(2,n2)),a=a,R=500)
  bci<-boot.ci(b,conf=conf,type="perc")
  Lci<-bci$percent[4]
  Uci<-bci$percent[5]
  return(list(LowerCI=Lci,UpperCI=Uci))
  
}

