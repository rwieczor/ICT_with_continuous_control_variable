# Example of ICT analysis


source("ICT_code.R")


# Example of survey data: 
# two columns y1 and y2 corresponding to X-Z and X+Z in our models
dat<-read.csv2("ankieta_full.csv")
y1<-dat$y1
y1<-na.omit(y1)
y2<-dat$y2
y2<-na.omit(y2)
n1<-length(y1)
n2<-length(y2)

# ML estimators:

em.algo.norm(y1,y2,pi_inits=NULL,mi_inits=NULL,sigma2_inits=NULL,
                   a=1,info=TRUE,AIC=TRUE)


em.algo.lognorm(y1,y2,pi_inits=NULL,mi_inits=NULL,sigma2_inits=NULL,
                   a=1,info=TRUE,AIC=TRUE)


em.algo.gamma(y1,y2,pi_inits=NULL,k_inits=NULL,theta_inits=NULL,
                   a=1,info=TRUE,AIC=TRUE)



# bootstrap confidence intervals:

ci.normal<-ML.boot.ci(y1,y2,a=1,conf=0.95,model="normal")
lapply(ci.normal,function(x) round(x,3))

ci.lognormal<-ML.boot.ci(y1,y2,a=1,conf=0.95,model="lognormal")
lapply(ci.lognormal,function(x) round(x,3))

ci.gamma<-ML.boot.ci(y1,y2,a=1,conf=0.95,model="gamma")
lapply(ci.gamma,function(x) round(x,3))



