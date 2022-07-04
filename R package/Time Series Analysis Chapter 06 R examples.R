setwd("C:/Wu Rongning/Teaching/CU/Statistics GR5221 & GU4221 Su2022/Lecture notes/Data sets") # Change path to where you save data file on your computer

#install.packages("tseries")
#install.packages("TSA")
#install.packages("sarima")
#install.packages("astsa")


# Section 6.1: ARIMA models for nonstationary time series

# --------------------------------------------------------------------
# Simulated ARIMA(1,1,0) process

set.seed(20)
arima.1=arima.sim(list(order=c(1,1,0),ar=0.8),n=200); arima.1=arima.1[2:201]
plot(arima.1,type="o",pch=22,lty=1,pty=2,xlab="",ylab="",main="Simulated ARIMA(1,1,0)")

op=par(mfrow=c(2,1))
acf(arima.1,lag.max=40,main="Simulated ARIMA(1,1,0)")
acf(arima.1,lag.max=40,type="partial",main="")
par(op)

# Lag-1 differencing
arima.1.diff=diff(arima.1)
plot(arima.1.diff,type="o",pch=22,lty=1,pty=2,xlab="",ylab="",main="Differenced ARIMA(1,1,0)")
abline(h=0)

op=par(mfrow=c(2,1))
acf(arima.1.diff,lag.max=40,main="Differenced ARIMA(1,1,0)")
acf(arima.1.diff,lag.max=40,type="partial",main="")
par(op)

# MLE AR(1) model for differenced series
ar(arima.1.diff,order.max=1,method="mle")

# MLE AR(2) model for simulated series
ar(arima.1,order.max=2,method="mle")

# MLE ARIMA(1,1,0) model for simulated series
arima(arima.1,order=c(1,1,0),method="ML")

# --------------------------------------------------------------------
# Dow-Jones Utilities Index, Aug 28-Dec 18, 1972; dowj.txt

dowj=read.table("dowj.txt",sep=""); dowj=dowj$V1
plot(dowj,type="o",pch=22,lty=1,pty=2,xlab="",ylab="",main="Dow-Jones Utilities Index, 8/28-12/18, 1972")

op=par(mfrow=c(2,1))
acf(dowj,lag.max=40,main="Dow-Jones Utilities Index")
acf(dowj,lag.max=40,type="partial",main="")
par(op)

# Lag-1 differencing
dowj.d=diff(dowj)
plot(dowj.d,type="o",pch=22,lty=1,pty=2,xlab="",main="Differenced Index")
abline(h=0)

op=par(mfrow=c(2,1))
acf(dowj.d,lag.max=40,main="Differenced Index")
acf(dowj.d,lag.max=40,type="partial",main="")
par(op)

# MLE AR(1) model for mean-corrected differenced data
ar(dowj.d,order.max=1,method="mle") # "demean=TRUE" by default

# MLE AR(2) model for original series
ar(dowj,order.max=2,method="mle")

# MLE ARIMA(1,1,0) model for original series
# Note: see 'Issue 2: your arima is drifting' in the file 'Some R time series issues.pdf'
arima(dowj,order=c(1,1,0),method="ML")
arima(dowj,order=c(1,1,0),xreg=1:length(dowj),method="ML")


# Section 6.2: Identification techniques

# --------------------------------------------------------------------
# Preliminary transformations: Australian red wine sales; wine.txt

wine=read.table("wine.txt",sep="")
wine=ts(wine,start=c(1980,1),frequency=12)
plot(wine,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="red wine sales")

# Log transformation
wine.log=log(wine)
plot(wine.log,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="log(red wine sales)")

# Differencing (at lag 12) log data
wine.log.d12=diff(wine.log,lag=12)
plot(wine.log.d12,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="differenced log(red wine sales)")
abline(h=0)

op=par(mfrow=c(2,1))
acf(wine.log.d12,lag.max=40,main="differenced log(red wine sales)")
acf(wine.log.d12,lag.max=40,type="partial",main="")
par(op)

# Fitting AR(12) model to the mean-corrected differenced data
arima(wine.log.d12-mean(wine.log.d12),order=c(12,0,0),include.mean=FALSE,method="ML")
ar.12=arima(wine.log.d12-mean(wine.log.d12),order=c(12,0,0),include.mean=FALSE,transform.pars=FALSE,fixed=c(NA,0,0,0,NA,0,0,0,0,0,0,NA),method="ML")
ar.12
# Note: see 'Issue 4: the wrong p-values' in the file 'Some R time series issues.pdf'

tsdiag(ar.12)


# Section 6.3: Unit roots in time series models

# --------------------------------------------------------------------
# Simulated trend stationary process & difference stationary process

# Trend stationary process
set.seed(10)
a0=5; a1=0.1; k=300
y1.sim=arima.sim(list(order=c(1,0,0),ar=0.8),n=k,rand.gen=rnorm,n.start=NA)
x1=a0+a1*seq(1,k)+y1.sim

# Difference stationary process
e.sim=rnorm(k,mean=0,sd=1)
y2.sim=diffinv(e.sim)
x2=a0+a1*seq(1,k)+y2.sim[2:301]

# Time series plots
x.min=min(x1,x2); x.max=max(x1,x2)
plot(seq(1,k),x1,type="l",pch=22,lty=2,pty=2,col="blue",ylim=c(x.min,x.max),xlab="",ylab="",main="")
lines(x2,col="red",lty=1)
legend(1,40,c("X_1", "X_2"),col=c("blue","red"),lty=c(2,1))

# Sample acf plots
op=par(mfrow=c(2,1))
acf(x1,lag.max=40,plot=TRUE,xlab="time lag",ylab="ACF",main="X_1")
acf(x2,lag.max=40,plot=TRUE,xlab="time lag",ylab="ACF",main="X_2")
par(op)

# --------------------------------------------------------------------
# Augmented Dickey-Fuller test; null hypothesis: x has an AR unit root

# Simulated ARIMA(1,1,0) process
set.seed(20)
arima.1=arima.sim(list(order=c(1,1,0),ar=0.8),n=200); arima.1=arima.1[2:201]

library(tseries)
adf.test(arima.1) # in "tseries" package
adf.test(arima.1,k=0) # standard Dickey-Fuller test

# Dow-Jones Utilities Index, Aug 28-Dec 18, 1972; dowj.txt
dowj=read.table("dowj.txt",sep=""); dowj=as.matrix(dowj)
adf.test(dowj)
adf.test(dowj,k=0)

# IID data
z=rnorm(1000)
adf.test(z) # no unit-root

# --------------------------------------------------------------------
# Phillips-Perron test; null hypothesis: x has an AR unit root

# Simulated ARIMA(1,1,0) process
library(tseries)
pp.test(arima.1) # in "tseries" package

# Dow-Jones Utilities Index, Aug 28-Dec 18, 1972; dowj.txt
pp.test(dowj)

# IID data
pp.test(z) # no unit-root

# --------------------------------------------------------------------
# Overdifferencing; data(rwalk) in R Package 'TSA'

library(TSA)
data(rwalk)
par(mfrow=c(2,1))
acf(diff(rwalk),lag.max=40,main="Differenced random walk")
acf(diff(rwalk,difference=2),lag.max=40,main="Twice differenced random walk")
detach("package:TSA", unload=TRUE)


# Section 6.4: Forecasting ARIMA models

# --------------------------------------------------------------------
# Dow-Jones Utilities Index, Aug 28-Dec 18, 1972; dowj.txt

dowj=read.table("dowj.txt",sep="")
dowj=as.ts(dowj$V1); (n.dowj=length(dowj))
plot(dowj,type="o",pch=22,lty=1,pty=2,xlab="",ylab="",main="Dow-Jones Utilities Index, 8/28-12/18, 1972")

op=par(mfrow=c(2,1))
acf(dowj,lag.max=40,main="Dow-Jones Utilities Index")
acf(dowj,lag.max=40,type="partial",main="")
par(op)

library(tseries)
adf.test(dowj)

# Fitting ARIMA(1,1,0) model to data using MLE
(fit.dowj=arima(dowj,order=c(1,1,0),xreg=1:n.dowj))
# Note: see 'Issue 2: your arima is drifting' in the file 'Some R time series issues.pdf'

# Diagnostic checking
tsdiag(fit.dowj)
# Note: see 'Issue 4: the wrong p-values' in the file 'Some R time series issues.pdf'

# Forecasting
n.pred=10
dowj.pred=predict(fit.dowj,n.pred,newxreg=(n.dowj+1):(n.dowj+n.pred))
ts.plot(dowj,dowj.pred$pred,col=1:2,xlab="",ylab="",main="Forecasting of Dow-Jones Utilities Index")
# Note: see 'Issue 3: you call that a forecast?' in the file 'Some R time series issues.pdf'


# Section 6.5: Seasonal ARIMA models

# --------------------------------------------------------------------
# ACF and PACF plots of SARMA models

set.seed(2021)
library(sarima)

# SMA models
sma.acf=ARMAacf(ar=,ma=c(rep(0,11),0.4),lag.max=50,pacf=FALSE)
sma.pacf=ARMAacf(ar=,ma=c(rep(0,11),0.4),lag.max=50,pacf=TRUE)
sma=sim_sarima(n=144,model=list(sma=0.4,nseasons=12,sigma2=1))
op=par(mfrow=c(1,3))
plot(seq(0,50),sma.acf,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="",ylim=c(-0.2,1.0))
abline(h=0)
plot(seq(0,50),rbind(c(1,sma.pacf)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="",ylim=c(-0.2,1.0))
abline(h=0)
plot(sma,type="o",pch=22,lty=1,pty=2,xlab="t",ylab="Simulated data")
par(op)

# SAR models
sar.acf=ARMAacf(ar=c(rep(0,11),0.5),ma=,lag.max=50,pacf=FALSE)
sar.pacf=ARMAacf(ar=c(rep(0,11),0.5),ma=,lag.max=50,pacf=TRUE)
sar=sim_sarima(n=144,model=list(sar=0.5,nseasons=12,sigma2=1))
op=par(mfrow=c(1,3))
plot(seq(0,50),sar.acf,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="",ylim=c(-0.2,1.0))
abline(h=0)
plot(seq(0,50),rbind(c(1,sar.pacf)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="",ylim=c(-0.2,1.0))
abline(h=0)
plot(sar,type="o",pch=22,lty=1,pty=2,xlab="t",ylab="Simulated data")
par(op)

# SARMA models
sarma1.acf=ARMAacf(ar=0.5,ma=c(rep(0,11),0.4),lag.max=50,pacf=FALSE)
sarma1.pacf=ARMAacf(ar=0.5,ma=c(rep(0,11),0.4),lag.max=50,pacf=TRUE)
sarma1=sim_sarima(n=144,model=list(ar=0.5,sma=0.4,nseasons=12,sigma2=1))
op=par(mfrow=c(1,3))
plot(seq(0,50),sarma1.acf,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="",ylim=c(-0.2,1.0))
abline(h=0)
plot(seq(0,50),rbind(c(1,sarma1.pacf)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="",ylim=c(-0.2,1.0))
abline(h=0)
plot(sarma1,type="o",pch=22,lty=1,pty=2,xlab="t",ylab="Simulated data")
par(op)

sarma2.acf=ARMAacf(ar=c(rep(0,11),0.5),ma=0.4,lag.max=50,pacf=FALSE)
sarma2.pacf=ARMAacf(ar=c(rep(0,11),0.5),ma=0.4,lag.max=50,pacf=TRUE)
sarma2=sim_sarima(n=144,model=list(sar=0.5,sma=0.4,nseasons=12,sigma2=1))
op=par(mfrow=c(1,3))
plot(seq(0,50),sarma2.acf,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="",ylim=c(-0.2,1.0))
abline(h=0)
plot(seq(0,50),rbind(c(1,sarma2.pacf)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="",ylim=c(-0.2,1.0))
abline(h=0)
plot(sarma2,type="o",pch=22,lty=1,pty=2,xlab="t",ylab="Simulated data")
par(op)

# --------------------------------------------------------------------
# Monthly accidental deaths data; deaths.txt

deaths=read.table("deaths.txt",sep=",")
deaths=ts(deaths,start=c(1973,1),frequency=12); (n=length(deaths))

# Differencing once at lag 12 and once lag 1; see plots below for justification
deaths.d12=diff(deaths,lag=12); deaths.d12.d1=diff(deaths.d12)

op=par(mfrow=c(2,2))
plot(deaths,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="monthly accidental deaths (in thousands)")
acf(deaths,lag.max=50,type="correlation",xlab="time lag",ylab="Sample ACF",ylim=c(-0.5,1),main="")
plot(deaths.d12,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="data after differencing at lag 12")
plot(deaths.d12.d1,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="data after differencings at lag 12 and lag 1")
par(op)

op=par(mfrow=c(2,1))
acf(deaths.d12.d1,lag.max=50,type="correlation",xlab="time lag",ylab="Sample ACF",ylim=c(-0.5,1),main="")
acf(deaths.d12.d1,lag.max=50,type="partial",xlab="time lag",ylab="Sample PACF",ylim=c(-0.5,1),main="")
par(op)

# Four candidate models for differenced data
(fit.1=arima(deaths.d12.d1-mean(deaths.d12.d1),order=c(0,0,1),seasonal=list(order=c(0,0,1),period=12),include.mean=FALSE))
arima(deaths.d12.d1-mean(deaths.d12.d1),order=c(1,0,0),seasonal=list(order=c(1,0,0),period=12),include.mean=FALSE)
arima(deaths.d12.d1-mean(deaths.d12.d1),order=c(0,0,1),seasonal=list(order=c(1,0,0),period=12),include.mean=FALSE)
arima(deaths.d12.d1-mean(deaths.d12.d1),order=c(1,0,0),seasonal=list(order=c(0,0,1),period=12),include.mean=FALSE)
tsdiag(fit.1)

# SARIMA model to original data and diagnostic checking
fit.deaths=arima(deaths,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12))
tsdiag(fit.deaths)

library(astsa)
sarima(deaths,0,1,1,0,1,1,12)

# Forecasting
n.pred=6
predict(fit.deaths,n.ahead=n.pred)

library(astsa)
sarima.for(deaths,n.ahead=n.pred,0,1,1,0,1,1,12)


# Section 6.6: Regression with ARMA errors

# --------------------------------------------------------------------
# Overshorts; oshorts.txt

oshorts=read.table("oshorts.txt",sep="")
overshorts=as.ts(oshorts$V1); (n.oshorts=length(overshorts))

plot(overshorts,type="o",pch=22,lty=1,pty=2,xlab="",ylab="(in gallons)",main="overshorts")
abline(h=0)

# ACF and PACF plots
op=par(mfrow=c(2,1))
acf(overshorts,lag.max=40,ylim=c(-0.5,1),main="overshorts")
acf(overshorts,lag.max=40,type="partial",ylim=c(-0.5,1),main="")
par(op)

# OLS estimate & its standard error, ignoring the dependence in the data
mean(overshorts); sd(overshorts)/sqrt(n.oshorts) # summary(glm(overshorts~1))

# MLE MA(1) model
arima(overshorts,order=c(0,0,1),method="ML") # 'include.mean=TRUE'

# --------------------------------------------------------------------
# Level of Lake Huron 1875-1972; lake.txt

lake=read.table("lake.txt",sep="")
lake=ts(lake$V1,start=1875); (n.lake=length(lake))

plot(lake,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="Lake Huron levels")

# ACF and PACF plots
op=par(mfrow=c(2,1))
acf(lake,lag.max=40,main="Lake Huron levels")
acf(lake,lag.max=40,type="partial",main="")
par(op)

# Least squares fit, ignoring the dependence in the data
trend=1:n.lake; summary(lm(lake~trend))

# MLE AR(2) model with linear trend
arima(lake,order=c(2,0,0),xreg=trend,method="ML")


