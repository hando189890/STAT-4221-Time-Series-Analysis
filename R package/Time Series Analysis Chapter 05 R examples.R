setwd("C:/Wu Rongning/Teaching/CU/Statistics GR5221 & GU4221 Su2022/Lecture notes/Data sets") # Change path to where you save data file on your computer

#install.packages("itsmr")
#install.packages("astsa")
#install.packages("tseries")
#install.packages("TSA")

# Section 5.1: Preliminary estimation

# --------------------------------------------------------------------
# Dow-Jones Utilities Index, Aug 28-Dec 18, 1972; dowj.txt

dowj=read.table("dowj.txt",sep=""); dowj=as.matrix(dowj)

op=par(mfrow=c(2,1))
plot(dowj,type="o",pch=22,lty=1,pty=2,xlab="",ylab="",main="Dow-Jones Utilities Index, 8/28-12/18, 1972")
acf(dowj,lag.max=40,main="Dow-Jones Utilities Index")
par(op)

# Lag-1 differencing
dowj.d=diff(dowj)
plot(dowj.d,type="o",pch=22,lty=1,pty=2,xlab="",main="Differenced Index")
abline(h=0)

op=par(mfrow=c(2,1))
acf(dowj.d,lag.max=40,main="Differenced Index")
acf(dowj.d,lag.max=40,type="partial",main="Differenced Index")
par(op)

mean(dowj.d)
dowj.d.mc=dowj.d-mean(dowj.d)

# Fitting AR(1) model to the mean-corrected differenced data using Yule-Walker
fit.dowj.1=ar(dowj.d.mc,order.max=1) # default method is "yule-walker"
fit.dowj.1$ar; fit.dowj.1$var.pred

# Fitting AR(1) model to the mean-corrected differenced data using Burg
ar(dowj.d.mc,order.max=1,method="burg")

# Fitting MA(2) model to the mean-corrected differenced data using innovations algorithm
library(itsmr)
ia(dowj.d.mc,2)

# --------------------------------------------------------------------
# El Nino and fish population; data(rec) in R package 'astsa'

library(astsa)
plot(rec,ylab="",xlab="",main="# of new fish")
acf2(rec-mean(rec),48) # producing values and graphic

rec.yw=ar.yw(rec-mean(rec),order=2,demean=FALSE)
rec.yw$ar # parameter estimates
sqrt(diag(rec.yw$asy.var.coef)) # standard errors of parameter estimates
rec.yw$var.pred # error variance estimate

acf(rec-mean(rec),2,type="covariance",plot=FALSE)$acf


# Section 5.2: Maximum likelihood estimation

# --------------------------------------------------------------------
# Level of Lake Huron 1875-1972; lake.txt

lake=read.table("lake.txt",sep="")
lake=ts(lake$V1,start=1875)
plot(lake,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="Lake Huron levels")

# Sample ACF and PACF plots
op=par(mfrow=c(2,1))
acf(lake,lag.max=40,main="Lake Huron levels")
acf(lake,lag.max=40,type="partial",main="")
par(op)

mean(lake)
x=lake-mean(lake)

# Fitting AR(2) to mean-corrected data and checking AIC
ar(x,method="mle",demean=FALSE)
ar(x,method="mle",demean=FALSE)$aic # AIC problem
arima(x,order=c(2,0,0),include.mean=FALSE,method="ML")

# Fitting ARMA(1,1) to mean-corrected data and checking AIC
arima(x,order=c(1,0,1),include.mean=FALSE,method="ML")

library(itsmr)
arma(x,p=1,q=1)

# Find the best model (with the lowest AICC) from a range of possible ARMA models
library(itsmr)
autofit(x,p=0:5,q=0:10) # slow when ranges of p and q are large


# Section 5.4: Diagnostic checking

# --------------------------------------------------------------------
# Level of Lake Huron 1875-1972; lake.txt

lake=read.table("lake.txt",sep="")
lake=ts(lake$V1,start=1875)

plot(lake,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="Lake Huron levels")

# Sample ACF and PACF plots
op=par(mfrow=c(2,1))
acf(lake,lag.max=40,main="Lake Huron levels")
acf(lake,lag.max=40,type="partial",main="")
par(op)

x=lake-mean(lake)

# Fitting ARMA(1,1) to mean-corrected data
arima(x,order=c(1,0,1),include.mean=FALSE,method="ML")

# Residuals from the fitted ARMA(1,1) model
y.res2=arima(x,order=c(1,0,1),include.mean=FALSE,method="ML")$residuals
plot(y.res2,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="Residuals from the fitted ARMA(1,1) model")
abline(h=0)

op=par(mfrow=c(2,1))
acf(y.res2,lag.max=40,main="Residuals from the fitted ARMA(1,1) model")
acf(y.res2,lag.max=40,type="partial",main="")
par(op)

Box.test(y.res2,lag=20,type="Box-Pierce")
Box.test(y.res2,lag=20,type="Ljung-Box")

library(TSA)
McLeod.Li.test(y=y.res2,gof.lag=20) # in package 'TSA'
detach("package:TSA", unload=TRUE)

qqnorm(y.res2); qqline(y.res2)
shapiro.test(y.res2)
library(tseries)
jarque.bera.test(y.res2) # in package 'tseries'


# Section 5.5: Forecasting

# --------------------------------------------------------------------
# Example: Overshorts; oshorts.txt

oshorts=read.table("oshorts.txt",sep="")
overshorts=as.ts(oshorts$V1); n=length(overshorts)

plot(overshorts,type="o",pch=22,lty=1,pty=2,xlab="",ylab="(in gallons)",main="overshorts")
abline(h=0)

op=par(mfrow=c(2,1))
acf(overshorts,lag.max=40,ylim=c(-0.5,1),main="overshorts")
acf(overshorts,lag.max=40,type="partial",ylim=c(-0.5,1),main="")
par(op)

mean(overshorts)

# Fitting MA(1) model to mean-corrected data
fit.oshorts=arima(overshorts-mean(overshorts),order=c(0,0,1),include.mean=FALSE,method="ML")
fit.oshorts$coef; fit.oshorts$sigma2

# Diagnostic checking
tsdiag(fit.oshorts)

# Forecasting
n.pred=5; (oshorts.pred=predict(fit.oshorts,n.pred))
ts.plot(overshorts,oshorts.pred$pred+mean(overshorts),col=1:2,xlab="",ylab="(in gallons)",main="Forecasting of overshorts")

# 95% prediction interval for X_62
mean(overshorts)+c(oshorts.pred$pred[n.pred]-1.96*oshorts.pred$se[n.pred],oshorts.pred$pred[n.pred]+1.96*oshorts.pred$se[n.pred])


# Section 5.6: Order selection

# --------------------------------------------------------------------
# Best subset ARMA selection based on BIC

set.seed(10)
test=arima.sim(model=list(ar=c(rep(0,11),0.5),ma=c(rep(0,11),0.7)),n=200)

library(TSA)
res=armasubsets(y=test,nar=14,nma=14,y.name='test',ar.method='ols')
detach("package:TSA", unload=TRUE)
par(mfrow=c(1,1)); plot(res)


