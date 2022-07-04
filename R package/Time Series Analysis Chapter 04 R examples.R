setwd("C:/Wu Rongning/Teaching/CU/Statistics GR5221 & GU4221 Su2022/Lecture notes/Data sets") # Change path to where you save data file on your computer

#install.packages("itsmr")

library(itsmr)

# --------------------------------------------------------------------
# Example: Sunspot numbers, 1770–1869; sunspots.txt

sunspots=read.table("sunspots.txt",sep="")
sun=ts(sunspots,start=1770)
plot(sun,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="count",main="sunspots")
mean(sun) # acvf.sun=acf(sun,lag.max=2,type="covariance",plot=FALSE); acvf.sun$acf[1:3]
abline(h=mean(sun))

op=par(mfrow=c(2,1))
acf(sun,lag.max=40,ylim=c(-0.7,1),main="Sample ACF of sunspots")
acf(sun,lag.max=40,type="partial",ylim=c(-0.7,1),main="Sample PACF of sunspots")
par(op)

ar(sun,aic=FALSE,order.max=2) # default method is "yule-walker"

plots(specify(ar=c(1.3175,-0.6341),sigma2=298.2)) # magnitude problem when applied to Problem 4.8(a)


# --------------------------------------------------------------------
# Example: spectral density of AR(1) models

set.seed(2021)

op=par(mfrow=c(2,3))
AR1.a=ARMAacf(ar=0.7,ma=,lag.max=20,pacf=FALSE)
plot(seq(0,20),AR1.a,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="AR(1), phi=0.7"); abline(h=0)
plots(specify(ar=0.7))
plot(arima.sim(n=200,list(ar=0.7)),ylab="a realization",main="AR(1), phi=0.7")

AR1.b=ARMAacf(ar=-0.7,ma=,lag.max=20,pacf=FALSE)
plot(seq(0,20),AR1.b,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="AR(1), phi=-0.7"); abline(h=0)
plots(specify(ar=-0.7))
plot(arima.sim(n=200,list(ar=-0.7)),ylab="a realization",main="AR(1), phi=-0.7")
par(op)


# --------------------------------------------------------------------
# Example: spectral density of MA(1) models

set.seed(2021)

op=par(mfrow=c(2,3))
MA1.a=ARMAacf(ar=,ma=0.9,lag.max=20,pacf=FALSE)
plot(seq(0,20),MA1.a,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="MA(1),theta=0.9"); abline(h=0)
plots(specify(ma=0.9))
plot(arima.sim(n=200,list(ma=0.9)),ylab="a realization",main="MA(1), theta=0.9")

MA1.b=ARMAacf(ar=,ma=-0.9,lag.max=20,pacf=FALSE)
plot(seq(0,20),MA1.b,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="MA(1), theta=-0.9"); abline(h=0)
plots(specify(ma=-0.9))
plot(arima.sim(n=200,list(ma=-0.9)),ylab="a realization",main="MA(1), theta=-0.9")
par(op)


# --------------------------------------------------------------------
# Example: periodogram of Sunspot numbers, 1770–1869; sunspots.txt

sunspots=read.table("sunspots.txt",sep="")
sun=ts(sunspots,start=1770)
plot(sun,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="count",main="sunspots"); abline(h=mean(sun))

op=par(mfrow=c(2,2))
plots(specify(ar=c(1.3175,-0.6341),sigma2=298.2))
periodogram(sun,q=0,opt=1) # Figure 4-9
periodogram(sun,q=1,opt=1) # Figure 4-10; 1 Daniell filter
periodogram(sun,q=c(1,2),opt=1) # Figure 4-11; 2 Daniell filters
par(op)


# --------------------------------------------------------------------
# Example: power transfer function of differencing filter

lambda=seq(0,pi,length=100)
ptf=2*(1-cos(12*lambda))
plot(lambda,ptf,type="l",main="")


