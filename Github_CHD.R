library(tikzDevice)

clr1<-rgb(27/255,158/255,119/255)
clr2<-rgb(217/255,95/255,2/255)
clr3<-rgb(117/255,112/255,179/255)
clr4<-rgb(231/255,41/255,138/255)
clr5<-rgb(102/255,166/255,30/255)
clr6<-rgb(230/255,171/255,2/255)
clr7<-rgb(166/255,118/255,29/255)


clr1b<-rgb(27/255,158/255,119/255,0.3)
clr2b<-rgb(217/255,95/255,2/255,0.3)
clr3b<-rgb(117/255,112/255,179/255,0.3)
clr4b<-rgb(231/255,41/255,138/255,0.3)
clr5b<-rgb(102/255,166/255,30/255,0.3)
clr6b<-rgb(230/255,171/255,2/255,0.3)
clr7b<-rgb(166/255,118/255,29/255,0.3)

clr1c<-rgb(27/255,158/255,119/255,0.5)
clr5c<-rgb(102/255,166/255,30/255,0.5)
clr6c<-rgb(230/255,171/255,2/255,0.5)


################################
#### Examples heterogeneity ####
################################


#BHN distributed U1

tikz("D:/Documents/PhD/Causal Hazard Difference/TikzFigures/FigBHN.tex",width=4,height=4)

plot(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){-log((0.5)*(exp(-x*1)-1)+1)}),col=clr2,type="l", ylim=c(-1,1),ylab="B(t)",xlab="t")
lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){-x*(-0.5)}),col=clr2b,lty=2)

lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){-log((0.5)*(exp(-x*(-0.25))-1)+(0.5)*(exp(-x*0)-1)+1)}),col=clr3)
lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){-x*(0.125)}),col=clr3b,lty=2)

lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){-log((0.5)*(exp(-x*(-0.1))-1)+(0.5)*(exp(-x*1)-1)+1)}),col=clr4)
#lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){-x*0.25}),col=clr3b,lty=2)
lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){-x*(-0.45)}),col=clr4b,lty=2)
abline(h=0,lwd=1,col=rgb(220/255,220/255,220/255,0.3))

legend(0.25, -0.4, title="(mu1, p1, mu2, p2)", legend=c("(0,0.5,1,0.5)","(-0.25,0.5,0,0.5)","(-0.1,0.5,1,0.5)"),lty=c(1,1,1),col=c(clr2,clr3,clr4),cex=0.75)
dev.off()

#Gamma distributed U1
tikz("D:/Documents/PhD/Causal Hazard Difference/TikzFigures/FigGamma.tex",width=4,height=4)

k<- 1
theta<-1
plot(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){k*log(theta*x+1)-0.5*x}),col=clr5,lwd=2,type="l",ylab="B(t)",xlab="t",ylim=c(-3,3))
lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){x*0.5}),col=clr5b,lty=1)


lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){k*log(theta*x+1)-1*x}),lwd=2,col=clr5,lty=2)
lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){x*0}),col=clr5b,lty=2)

lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){k*log(theta*x+1)-0.25*x}),lwd=2,col=clr5,lty=3)
lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){x*0.75}),col=clr5b,lty=3)


lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){k*log(theta*x+1)}),lwd=2,col=clr5,lty=4)
lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){x*1}),col=clr5b,lty=4)


legend(0.5, -1, title="k", legend=c(0, 0.25, 0.5, 1),lty=c(4,3,1,2),col=c(clr5,clr5,clr5,clr5),seg.len=3, lwd=2,cex=0.75)
dev.off()

##################
####Dependence####
##################

##B(t) curve

library(copula)
library(gumbel)
library(GoFKernel)
library(Compounding)
library(stabledist)

g <- function(a) {
    iu <- complex(real=0, imaginary=1)
    return(abs(1 - iu * tan(pi * a / 2)) ^ (-1 / a))
}

library(statmod)


simcop<-function(n, cop, tau, frail, var0=1, alpha=0.5, HR, var1=1, k){
    
    if(cop=="Gaussian"){
        #####Gaussian########
        rho<- sin(pi*tau/2)
        nc<- normalCopula(rho)
        NU  <- rCopula(n, copula = nc)
    }
    
    
    #####Frank Copula####
    if(cop=="Frank"){
        theta <- iTau(copFrank, tau) # copula parameter
        fc <- frankCopula(theta, dim = 2) # define a Frank copula
        NU  <- rCopula(n, copula = fc)
    }
    
    
    #####Gumbel #####
    if(cop=="Gumbel"){
        theta<-iTau(copGumbel, tau)
        
        if(theta>=1){
            gc <- gumbelCopula(theta) # (note the default dim = 2)
            NU  <- rCopula(n, copula = gc)
        }
        
        if(theta<1){
            theta<-iTau(copGumbel, -tau)
            gc <- gumbelCopula(theta) # (note the default dim = 2)
            NU  <- rCopula(n, copula = gc)
            NU[,1]<-1-NU[,1]
        }
        
    }
    
    
    ####Clayton ######
    if(cop=="Clayton"){
        theta<-iTau(copClayton, tau)
        cc <- claytonCopula(theta, dim = 2) 
        
        NU  <- rCopula(n, copula = cc)
    }
    
    if(cop=="Independent"){
        NU  <- cbind(runif(n),runif(n))
    }
    
    if(cop=="Dependent"){
        NU  <- runif(n)
        if(tau==1){
            NU <-cbind(NU,NU)
        }
        if(tau== -1){
            NU <-cbind(NU,1-NU)
        }
    }
    
    ########################
    ###Simulate U0, U1 #####
    ########################
    
    if(frail=="gamma"){
        U0<-sapply(NU[,1],function(x){qgamma(x, shape=1/var0,scale=var0)})
        
    }
    
    if(frail=="cpois"){
        cpoisdata<-(sapply(rpois(n,3/var0),function(x){ifelse(x>0,sum(rgamma(x,shape=1/2,scale=(2/3)*var0)),0)}))
        U0<-quantile(cpoisdata,NU[,1])
        
    }
    
    if(frail=="pstable"){
        gamma<-g(alpha)
        U0<-sapply(NU[,1],function(x){tryCatch({qstable(x, alpha=alpha, beta=1, gamma=gamma, delta=0, pm=1)}, error=function(e){qstable(round(x,3), alpha=alpha, beta=1, gamma=gamma, delta=0, pm=1)})})
        
    }
    
    if(frail=="invgauss"){
       
        
        U0<-sapply(NU[,1],function(x){qinvgauss(x,1,1/var0)})    
    }
    
   
    
    U1<-sapply(NU[,2],function(x){qgamma(x, shape=HR^2/var1,scale=var1/HR)})-k
    
    NT<-runif(n,0,1)
    
    T0<-(9*(U0/3.000)^2*log(1/NT)+sqrt(3)*sqrt(27*(U0/3.000)^4*log(1/NT)^2+4*(U0/3.000)^3*(k+0)^3))^(1/3)/(2^(1/3)*3^(2/3)*(U0/3.000))-((k+0)*(2/3)^(1/3))/((9*(U0/3.000)^2*log(1/NT)+sqrt(3)*sqrt(27*(U0/3.000)^4*log(1/NT)^2+4*(U0/3.000)^3*(k+0)^3))^(1/3))
    T1<-(9*(U0/3.000)^2*log(1/NT)+sqrt(3)*sqrt(27*(U0/3.000)^4*log(1/NT)^2+4*(U0/3.000)^3*(k+U1)^3))^(1/3)/(2^(1/3)*3^(2/3)*(U0/3.000))-((k+U1)*(2/3)^(1/3))/((9*(U0/3.000)^2*log(1/NT)+sqrt(3)*sqrt(27*(U0/3.000)^4*log(1/NT)^2+4*(U0/3.000)^3*(k+U1)^3))^(1/3))

    # T0<-(9*(U0/60.000)^2*log(1/NT)+sqrt(3)*sqrt(27*(U0/60.000)^4*log(1/NT)^2+4*(U0/60.000)^3*(k+0)^3))^(1/3)/(2^(1/3)*3^(2/3)*(U0/60.000))-((k+0)*(2/3)^(1/3))/((9*(U0/60.000)^2*log(1/NT)+sqrt(3)*sqrt(27*(U0/60.000)^4*log(1/NT)^2+4*(U0/60.000)^3*(k+0)^3))^(1/3))
    # T1<-(9*(U0/60.000)^2*log(1/NT)+sqrt(3)*sqrt(27*(U0/60.000)^4*log(1/NT)^2+4*(U0/60.000)^3*(k+U1)^3))^(1/3)/(2^(1/3)*3^(2/3)*(U0/60.000))-((k+U1)*(2/3)^(1/3))/((9*(U0/60.000)^2*log(1/NT)+sqrt(3)*sqrt(27*(U0/60.000)^4*log(1/NT)^2+4*(U0/60.000)^3*(k+U1)^3))^(1/3))
    
 
    pA<-0.5
    A<-rbinom(n,1,pA)
    
    T<-A*T1+(1-A)*T0
    
    return(data.frame(A=A, T=T, T0=T0, T1=T1, U0=U0, U1=U1))
    
}

#cop: "Gaussian", "Frank", "Gumbel", "Clayton"
#frail: "gamma", "cpois", "pstable", "invgauss"
par(mfrow = c(1,3))

n<-10000000
cop<-"Gaussian"
frail<-"gamma"

tikz("D:/Documents/PhD/Causal Hazard Difference/TikzFigures/FigCopNew.tex",width=12,height=4)
par(mfrow = c(1,3))

for(k in c(0,0.5,1)){
    var0<-1
    var1<-1
    
    
    
    dataH0<-simcop(n, cop='Independent', tau= -1, frail=frail, var0=var0, alpha=0, HR=1, var1=var1, k)
    dataHneg<-simcop(n, cop='Dependent', tau= -1, frail=frail, var0=var0, alpha=0, HR=1, var1=var1, k)
    dataHpos<-simcop(n, cop='Dependent', tau= 1, frail=frail, var0=var0, alpha=0, HR=1, var1=var1, k)
    
    U1<-dataH0$U1
    U0<-dataH0$U0
    T1<-dataH0$T1
    T0<-dataH0$T0  
    
    A<-dataH0$A
    T<-dataH0$T
    
    #Other scenario: change 1.000 to 20.000
    HD<-sapply(seq(0.1,10,0.1),function(t){mean(U1[which(T1>t)])+mean(U0[which(T1>t)]-U0[which(T1>t)])*(t^2/1.000)})

    if(k==0){
    plot(seq(0.1,10,0.1),cumsum(HD)*0.1,ylim=c(0,5),xlim=c(0,5),ylab="B(t)", xlab="t",col=clr5c,type="l",lty=1,axes=T,cex.lab=1.5,cex.axis=1.5)
    }
    if(k==0.5){
        plot(seq(0.1,10,0.1),cumsum(HD)*0.1,ylim=c(-2,3),xlim=c(0,5),ylab="B(t)", xlab="t",col=clr5c,type="l",lty=1,axes=T,cex.lab=1.5,cex.axis=1.5)
    }
    if(k==1){
    plot(seq(0.1,10,0.1),cumsum(HD)*0.1,ylim=c(-4.5,0.5),xlim=c(0,5),ylab="B(t)", xlab="t",col=clr5c,type="l",lty=1,axes=T,cex.lab=1.5,cex.axis=1.5)
        
    }
    
    #For survival curves
    #plot(seq(0,5,0.1),sapply(seq(0,5,0.1),function(t){length(which(T0>t))/n}),ylim=c(0,1),ylab="S(t)", xlab="t",col="black",type="l",lty=1,axes=T,cex.lab=1.5,cex.axis=1.5)
     
     
     
    lines(seq(0,5,0.1),sapply(seq(0,5,0.1),function(x){x*(1-k)}),col="lightgrey",lty=1)

    
    U1<-dataHneg$U1
    U0<-dataHneg$U0
    T1<-dataHneg$T1
    T0<-dataHneg$T0
    
    #Other scenario: change 1.000 to 20.000
    HD<-sapply(seq(0.1,10,0.1),function(t){mean(U1[which(T1>t)])+mean(U0[which(T1>t)]-U0[which(T1>t)])*(t^2/1.000)})
    lines(seq(0.1,10,0.1),cumsum(HD)*0.1,ylab="Hazard ratio",xlab="t",col=clr6c,ylim=c(0,4),lty=5)
    
    type<-4
    for(tau in c(-0.5, 0.5)){
        data<-simcop(n, cop=cop, tau= tau, frail=frail, var0=var0, alpha=0, HR=1, var1=var1, k)
        
        U1<-data$U1
        U0<-data$U0
        T1<-data$T1
        T0<-data$T0
        
        #Other scenario: change 1.000 to 20.000
        HD<-sapply(seq(0.1,10,0.1),function(t){mean(U1[which(T1>t)])+mean(U0[which(T1>t)]-U0[which(T1>t)])*(t^2/1.000)})
        lines(seq(0.1,10,0.1),cumsum(HD)*0.1,ylab="Hazard ratio",xlab="t",col=clr6c,ylim=c(0,4),lty=type)

        type<-3
    }
    
    
    U1<-dataHpos$U1
    U0<-dataHpos$U0
    T1<-dataHpos$T1
    T0<-dataHpos$T0
    
    #Other scenario: change 1.000 to 20.000
    HD<-sapply(seq(0.1,10,0.1),function(t){mean(U1[which(T1>t)])+mean(U0[which(T1>t)]-U0[which(T1>t)])*(t^2/1.000)})
    lines(seq(0.1,10,0.1),cumsum(HD)*0.1,ylab="Hazard ratio",xlab="t",col=clr6c,ylim=c(0,4),lty=2)

    U1<-dataH0$U1
    U0<-dataH0$U0
    T1<-dataH0$T1
    T0<-dataH0$T0   
    
    #Other scenario: change 1.000 to 20.000
    HD<-sapply(seq(0.1,10,0.1),function(t){mean(U1[which(T1>t)])+mean(U0[which(T1>t)]-U0[which(T1>t)])*(t^2/1.000)})
    lines(seq(0.1,10,0.1),cumsum(HD)*0.1,ylim=c(-3,3),xlim=c(0,5),ylab="Hazard ratio", xlab="t",col=clr5c,type="l",lty=1)
    
    if(k==0){
        legend(0.5, 4.5, title="rho", legend=c(-1, 0.5, 0, 0.5, 1),lty=c(5,4,1,3,2),col=c(clr6c,clr6c,clr5c,clr6c,clr6c),seg.len=3, lwd=2,cex=1.5)
    }
    
}
dev.off()


################
## Case study ##
################

library("invGauss")
data(d.oropha.rec)

library(timereg)
aalenfit<-aalen(Surv(time,status)~factor(treatm),data=d.oropha.rec[d.oropha.rec$site%in%c(1,4),])

tikz("D:/Documents/PhD/Causal Hazard Difference/TikzFigures/CaseStudy1b.tex",width=4,height=4)
plot(aalenfit,ylim=c(-1,1),main=F,lwd=1)

p1<- 0.5
mu1<- -0.1
mu2<- 0.4

lines(seq(0,4.3,0.01),sapply(seq(0,4.3,0.01),function(x){-log((p1)*(exp(-x*mu1)-1)+(1-p1)*(exp(-x*mu2)-1)+1)}),col=clr1,lwd=2)
dev.off()

##Simulation##
# n<-130
# 
# p1<-0.5
# mu1<- -0.1
# mu2<- 0.4
# 
# U1p<-rbinom(n,1,p1)
# U1<-mu1*U1p+mu2*(1-U1p)
# NT<-runif(n,0,1)
# T0<-log(1/NT)/(0.45)
# T1<-log(1/NT)/(0.45+U1)
# 
# 
# pA<-0.5
# A<-rbinom(n,1,pA)
# 
# T<-A*T1+(1-A)*T0
# 
# test.data<-data.frame(A=A, T=T, event=1)
# 
# 
# test.aalenfit<-aalen(Surv(T,event)~factor(A),data=test.data)
# 
# plot(test.aalenfit,xlim=c(0,4),ylim=c(-1,1))
# lines(seq(0,10,0.01),sapply(seq(0,10,0.01),function(x){-log((p1)*(exp(-x*mu1)-1)+(1-p1)*(exp(-x*mu2)-1)+1)}),col=clr1)

#Modification by site of tumor
tikz("D:/Documents/PhD/Causal Hazard Difference/TikzFigures/CaseStudy2.tex",width=8,height=4)
#plot(aalenfit,ylim=c(-1,1),main=F,lwd=1)
par(mfrow=c(1,2))
aalenfit2a<-aalen(Surv(time,status)~factor(treatm),data=d.oropha.rec[d.oropha.rec$site==1,])
plot(aalenfit2a,specific.comps=2,main=F)
lines(seq(0,4,0.01),sapply(seq(0,4,0.01),function(x){x*(0.4)}),col=clr1c)

aalenfit2c<-aalen(Surv(time,status)~factor(treatm),data=d.oropha.rec[d.oropha.rec$site==4,])
plot(aalenfit2c,specific.comps=2,main=F)
lines(seq(0,4,0.01),sapply(seq(0,4,0.01),function(x){x*(-0.1)}),col=clr1c)
dev.off()




