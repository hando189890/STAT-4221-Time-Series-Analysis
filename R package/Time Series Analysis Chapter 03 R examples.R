setwd("C:/Wu Rongning/Teaching/CU/Statistics GR5221 & GU4221 Su2022/Lecture notes/Data sets") # Change path to where you save data file on your computer

# --------------------------------------------------------------------
# Convert ARMA process to infinite MA process

ARMAtoMA(ar=c(0.7,-0.1),ma=,lag.max=20) # AR(2) process X_{t}-0.7X_{t-1}+0.1X_{t-2}=Z_{t}

ARMAtoMA(ar=0.5,ma=0.4,lag.max=10) # ARMA(1,1) process X_{t}-0.5X_{t-1}=Z_{t}+0.4Z_{t-1}


# --------------------------------------------------------------------
# ACF of some MA(q) processes

h=20

op=par(mfrow=c(2,3))
acf.MA1=ARMAacf(ar=,ma=0.5,lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.MA1,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(a) Xt = (1+0.5B)Zt",ylim=c(-0.2,1.0))
abline(h=0)

acf.MA2=ARMAacf(ar=,ma=c(0.3,-0.1),lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.MA2,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(b) Xt = (1+0.3B-0.1B^2)Zt",ylim=c(-0.2,1.0))
abline(h=0)

acf.MA3.1=ARMAacf(ar=,ma=c(1,0.11,-0.07),lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.MA3.1,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(c1) Xt = (1+B+0.11B^2-0.07B^3)Zt",ylim=c(-0.2,1.0))
abline(h=0)

acf.MA3.2=ARMAacf(ar=,ma=c(1,0.31,0.3),lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.MA3.2,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(c2) Xt = (1+B+0.31B^2+0.3B^3)Zt",ylim=c(-0.2,1.0))
abline(h=0)

acf.MA4.1=ARMAacf(ar=,ma=c(0.7,-0.19,-0.103,0.021),lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.MA4.1,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(d1) Xt = (1+0.7B-0.19B^2-0.103B^3+0.021B^4)Zt",ylim=c(-0.2,1.0))
abline(h=0)

acf.MA4.2=ARMAacf(ar=,ma=c(1.1,0.41,0.061,0.003),lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.MA4.2,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(d2) Xt = (1+1.1B+0.41B^2+0.061B^3+0.003B^4)Zt",ylim=c(-0.2,1.0))
abline(h=0)
par(op)


# --------------------------------------------------------------------
# Overshorts; oshorts.txt

oshorts=read.table("oshorts.txt",sep="")
oshorts=as.matrix(oshorts)

op=par(mfrow=c(2,1))
plot(oshorts,type="o",pch=22,lty=1,pty=2,xlab="",ylab="(in gallons)",main="overshorts"); abline(h=0)
acf(oshorts,lag.max=40,main="Sample ACF of overshorts")
par(op)


# --------------------------------------------------------------------
# ACF of some AR(2) processes

h=20

op=par(mfrow=c(2,2))
acf.AR2.a=ARMAacf(ar=c(0.7,-0.1),ma=,lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.AR2.a,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(a) (1-0.5B)(1-0.2B)Xt = Zt",ylim=c(-1.0,1.0))
abline(h=0)

acf.AR2.b=ARMAacf(ar=c(1.4,-0.45),ma=,lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.AR2.b,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(b) (1-0.9B)(1-0.5B)Xt = Zt",ylim=c(-1.0,1.0))
abline(h=0)

acf.AR2.c=ARMAacf(ar=c(-0.4,0.45),ma=,lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.AR2.c,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(c) (1+0.9B)(1-0.5B)Xt = Zt",ylim=c(-1.0,1.0))
abline(h=0)

acf.AR2.d=ARMAacf(ar=c(0.75,-0.5625),ma=,lag.max=h,pacf=FALSE)
plot(seq(0,h),acf.AR2.d,type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="ACF",main="(d) (1-0.75B+0.5625B^2)Xt = Zt",ylim=c(-1.0,1.0))
abline(h=0)
par(op)


# --------------------------------------------------------------------
# PACF of some AR(p) processes

h=20

op=par(mfrow=c(2,2))
pacf.AR1=ARMAacf(ar=0.5,ma=,lag.max=h,pacf=TRUE)
plot(seq(0,h),rbind(c(1,pacf.AR1)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="(a) (1-0.5B)Xt=Zt",ylim=c(-0.5,1.0))
abline(h=0)

pacf.AR2=ARMAacf(ar=c(0.7,-0.1),ma=,lag.max=h,pacf=TRUE)
plot(seq(0,h),rbind(c(1,pacf.AR2)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="(b) (1-0.5B)(1-0.2B)Xt=Zt",ylim=c(-0.5,1.0))
abline(h=0)

pacf.AR3=ARMAacf(ar=c(0.5,0.11,-0.07),ma=,lag.max=h,pacf=TRUE)
plot(seq(0,h),rbind(c(1,pacf.AR3)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="(c) (1-0.5B-0.11B^2+0.07B^3)Xt=Zt",ylim=c(-0.5,1.0))
abline(h=0)

pacf.AR4=ARMAacf(ar=c(0.7,-0.19,-0.103,0.021),ma=,lag.max=h,pacf=TRUE)
plot(seq(0,h),rbind(c(1,pacf.AR4)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="(d) (1-0.7B+0.19B^2+0.103B^3-0.021B^4)Xt=Zt",ylim=c(-0.5,1.0))
abline(h=0)
par(op)


# --------------------------------------------------------------------
# PACF of some MA(q) processes

h=20

op=par(mfrow=c(2,2))
pacf.MA1=ARMAacf(ar=,ma=0.5,lag.max=h,pacf=TRUE)
plot(seq(0,h),rbind(c(1,pacf.MA1)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="(a) Xt=(1+0.5B)Zt",ylim=c(-0.5,1.0))
abline(h=0)

pacf.MA2=ARMAacf(ar=,ma=c(0.3,-0.1),lag.max=h,pacf=TRUE)
plot(seq(0,h),rbind(c(1,pacf.MA2)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="(b) Xt=(1+0.3B-0.1B^2)Zt",ylim=c(-0.5,1.0))
abline(h=0)

pacf.MA3=ARMAacf(ar=,ma=c(1,0.11,-0.07),lag.max=h,pacf=TRUE)
plot(seq(0,h),rbind(c(1,pacf.MA3)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="(c) Xt=(1+B+0.11B^2-0.07B^3)Zt",ylim=c(-0.5,1.0))
abline(h=0)

pacf.MA4=ARMAacf(ar=,ma=c(0.7,-0.19,-0.103,0.021),lag.max=h,pacf=TRUE)
plot(seq(0,h),rbind(c(1,pacf.MA4)),type="h",pch=22,lty=1,pty=2,xlab="Lag",ylab="PACF",main="(d) Xt=(1+0.7B-0.19B^2-0.103B^3+0.021B^4)Zt",ylim=c(-0.5,1.0))
abline(h=0)
par(op)


# --------------------------------------------------------------------
# Sunspot numbers; sunspots.txt

sunspots=read.table("sunspots.txt",sep="")
sunspots=as.matrix(sunspots); n=length(sunspots)

year=seq(1770,by=1,length=n)
op=par(mfrow=c(2,1))
plot(year,sunspots,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="count",main="sunspots"); abline(h=mean(sunspots))
acf(sunspots,lag.max=40,type="partial",main="Sample PACF of sunspots")
par(op)


