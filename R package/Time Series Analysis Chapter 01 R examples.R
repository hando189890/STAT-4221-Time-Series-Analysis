setwd("C:/Wu Rongning/Teaching/CU/Statistics GR5221 & GU4221 Su2022/Lecture notes/Data sets") # Change path to where you save data file on your computer

#install.packages("astsa")
#install.packages("TSA")
#install.packages("tseries")


# Section 1.1: Examples of time series

# --------------------------------------------------------------------
# Australian red wine sales; wine.txt

wine=read.table("wine.txt",sep=",")
wine=ts(wine,start=c(1980,1),frequency=12)
plot(wine,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="red wine sales (in kiloliters)")
# Explanation of pch can be found in ?points

# --------------------------------------------------------------------
# Accidental deaths, USA, 1973-1978; deaths.txt

deaths=read.table("deaths.txt",sep=",")
deaths=ts(deaths,start=c(1973,1),frequency=12)
plot(deaths,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="monthly accidental deaths (in thousands)")

# --------------------------------------------------------------------
# Dow-Jones Index (closing prices; 9/10/93-8/26/94); dowj2.csv

dj2=read.csv("dowj2.csv",header=FALSE,sep=",")
dj2=as.matrix(dj2)
plot(dj2,type="o",pch=22,lty=1,pty=2,xlab="",ylab="Closing prices")

# --------------------------------------------------------------------
# Population of the USA, 1790-1990; uspop.txt

uspop=read.table("uspop.txt",sep=",")
uspop=as.matrix(uspop); n=length(uspop)

year=seq(1790,by=10,length=n); #length(year)
plot(year,uspop/10^6,type="o",pch=22,lty=1,pty=2,xaxt="n",xlab="year",ylab="population of the USA (in millions)")
# Explanation of xaxt can be found in ?par
axis(1,at=seq(1800,1980,by=20),labels=c("1800","1820","1840","1860","1880","1900","1920","1940","1960","1980"))


# Section 1.3: Simple time series models

# --------------------------------------------------------------------
# Population of the USA, 1790-1990; uspop.txt

uspop=read.table("uspop.txt",sep=",")
uspop=as.matrix(uspop); n=length(uspop)

year=seq(1790,by=10,length=n)
plot(year,uspop/10^6,type="o",pch=22,lty=1,pty=2,xaxt="n",xlab="year",ylab="population of the USA (in millions)")
axis(1,at=seq(1800,1980,by=20),labels=c("1800","1820","1840","1860","1880","1900","1920","1940","1960","1980"))

# Least squares fit
year.2=year^2
data.uspop=data.frame(uspop,year,year.2)
fit.uspop=glm(uspop~1+year+year.2,data=data.uspop)
fit.uspop$coefficients; #summary(fit.uspop); fit.uspop$fitted.values
lines(year,fit.uspop$fitted.values/10^6,lty=1,col="red")

# --------------------------------------------------------------------
# Accidental deaths, USA, 1973-1978; deaths.txt

deaths=read.table("deaths.txt",sep="")
deaths=as.matrix(deaths); n.d=length(deaths)

year=seq(1973,by=1/12,length=n.d)
plot(year,deaths/1000,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="monthly accidental deaths (in thousands)")

# Harmonic regression
t=seq(1,n.d)
lambda.1=2*pi/12; lambda.2=2*pi/6
cos.1=cos(lambda.1*t); cos.2=cos(lambda.2*t); sin.1=sin(lambda.1*t); sin.2=sin(lambda.2*t)
data.deaths=data.frame(deaths,cos.1,cos.2,sin.1,sin.2)

fit.deaths=lm(deaths~1+cos.1+cos.2+sin.1+sin.2,data=data.deaths)
summary(fit.deaths)

plot(year,deaths/1000,type="p",pch=22,lty=1,pty=2,xlab="year",ylab="monthly accidental deaths (in thousands)")
lines(year,fit.deaths$fitted.values/1000,lty=1,col="red")


# Section 1.4: Stationary models and ACF

# --------------------------------------------------------------------
# Sample ACF of 200 simulated values of iid N(0,1) noise

set.seed(49)
m=200; time=seq(1,200,by=1)
x=rnorm(m,mean=0,sd=1)
plot(time,x,type="o",pch=22,lty=1,pty=2,xlab="",ylab="")
abline(h=0)

acf(x,lag.max=40,type="correlation",plot=TRUE,xlab="time lag",ylab="Sample ACF",ylim=c(-0.2,1),main="")

# --------------------------------------------------------------------
# Sample ACF of wine.txt

wine=read.table("wine.txt",sep=",")
wine=as.matrix(wine); n=length(wine)

year=seq(1980,by=1/12,length=n)

par(mfrow=c(1,2))
plot(year,wine,type="o",pch=22,lty=1,pty=2,	xlab="year",ylab="red wine sales (in kiloliters)")
acf(wine,lag.max=40,type="correlation",plot=TRUE,xlab="time lag",ylab="ACF",ylim=c(-0.2,1),main="")


# Section 1.5: Trend and seasonality removal

# --------------------------------------------------------------------
# Box-Cox transformation: Australian red wine sales; wine.txt

wine=read.table("wine.txt",sep="")
wine=ts(wine,start=c(1980,1),frequency=12)
plot(wine,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="red wine sales")

# Log transformation
wine.log=log(wine)
plot(wine.log,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="log(red wine sales)")

# --------------------------------------------------------------------
# Classical decomposition: monthly accidental deaths data; deaths.txt

# Additive decomposition 
deaths=read.table("deaths.txt",sep="")
deaths.ts=ts(deaths,frequency=12,start=c(1973,1))
plot(deaths.ts/1000,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="monthly accidental deaths (in thousands)")

deaths.de=decompose(deaths.ts,type="additive")

op=par(mfrow=c(3,1))
plot(deaths.de$seasonal/1000,type="b",pch=22,lty=1,pty=2,xlab="year",ylab="seasonal component")
plot(deaths.de$trend/1000,type="b",pch=22,lty=1,pty=2,xlab="year",ylab="trend component")
plot(deaths.de$random/1000,type="b",pch=22,lty=1,pty=2,xlab="year",ylab="random component")
par(op)

# --------------------------------------------------------------------
# Weekly mortality in Los Angeles County; data(cmort) in R package 'astsa'

library(astsa)
data(package='astsa')

# Polynomial and periodic regression smoothers
wk=time(cmort)-mean(time(cmort)) # wk is essentially t/52 centered at zero
wk2=wk^2; wk3=wk^3
cs=cos(2*pi*wk); sn=sin(2*pi*wk)
reg1=lm(cmort~wk+wk2+wk3,na.action=NULL)
reg2=lm(cmort~wk+wk2+wk3+cs+sn,na.action=NULL)
plot(cmort,type="p",ylab="mortality")
lines(fitted(reg1),col="red"); lines(fitted(reg2),col="blue")

# Moving average smoother
ma5=filter(cmort,sides=2,rep(1,5)/5)
ma53=filter(cmort,sides=2,rep(1,53)/53)
plot(cmort,type="p",ylab="mortality")
lines(ma5,col="blue"); lines(ma53,col="red")

# Kernel smoother
plot(cmort,type="p",ylab="mortality")
lines(ksmooth(time(cmort),cmort,"normal",bandwidth=5/52),col="blue")
lines(ksmooth(time(cmort),cmort,"normal",bandwidth=2),col="red")

# --------------------------------------------------------------------
# Differencing: monthly accidental deaths data; deaths.txt

deaths=read.table("deaths.txt",sep="")
deaths.ts=ts(deaths,frequency=12,start=c(1973,1))
plot(deaths.ts/1000,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="monthly accidental deaths (in thousands)")

# Deseasonalized monthly accidental deaths
deaths.d=diff(deaths.ts,lag=12)
plot(deaths.d/1000,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="(thousands)")

# Detrended and deseasonalized monthly accidental deaths
deaths.dd=diff(deaths.d)
plot(deaths.dd/1000,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="(thousands)")
abline(h=0)


# Section 1.6: Testing the estimated noise sequence

# --------------------------------------------------------------------
# 200 simulated values of iid N(0,1) noise

set.seed(49)
m=200; time=seq(1,200,by=1)
x=rnorm(m,mean=0,sd=1)
plot(time,x,type="o",pch=22,lty=1,pty=2,xlab="",ylab="")
abline(h=0)

acf(x,lag.max=40,xlab="time lag",ylab="Sample ACF",ylim=c(-0.2,1),main="")

Box.test(x,lag=20,type="Box-Pierce")
Box.test(x,lag=20,type="Ljung-Box")

library(TSA)
McLeod.Li.test(y=x,gof.lag=20) # in package 'TSA'
detach("package:TSA",unload=TRUE)

qqnorm(x); qqline(x)
shapiro.test(x)

library(tseries)
jarque.bera.test(x) # in package 'tseries'

# --------------------------------------------------------------------
# Level of Lake Huron 1875-1972; lake.txt

lake=read.table("lake.txt",sep="")
lake=ts(lake,start=1875)

plot(lake,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="Lake Huron levels")

# Least squares fit
fit.lake=lm(lake~1+time(lake))
lines((1875:1972),fit.lake$fitted.values,lty=1,col="red")

# Residual plot
y=fit.lake$residuals
plot(y,type="o",pch=22,lty=1,pty=2,xlab="year",ylab="",main="Residuals from a linear fit")
abline(h=0)

# ACF plot of residual series
acf(y,lag.max=40,main="Lake Huron residuals from a linear fit")

Box.test(y,lag=20,type="Box-Pierce")
Box.test(y,lag=20,type="Ljung-Box")

library(TSA)
McLeod.Li.test(y=y,gof.lag=20) # in package 'TSA'
detach("package:TSA",unload=TRUE)

qqnorm(y); qqline(y)
shapiro.test(y)
jarque.bera.test(y) # in package 'tseries'


