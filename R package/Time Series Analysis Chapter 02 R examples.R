setwd("C:/Wu Rongning/Teaching/CU/Statistics GR5221 & GU4221 Su2022/Lecture notes/Data sets") # Change path to where you save data file on your computer

#install.packages("TSA")

# --------------------------------------------------------------------
# A simulated MA(1) series with coefficient equal to -0.9; data(ma1.1.s) in R package 'TSA'

library(TSA)
data(ma1.1.s)

plot(ma1.1.s,xlab="",ylab="",main="A simulated MA(1) series"); abline(h=0)

par(mfrow=c(2,1))
acf(ma1.1.s,lag.max=40,main="A simulated MA(1) series")
acf(ma1.1.s,ci.type='ma',lag.max=40,main="A simulated MA(1) series")

detach("package:TSA",unload=TRUE)


# --------------------------------------------------------------------
# Level of Lake Huron 1875-1972; lake.txt

huron=read.table("lake.txt",sep=",")
huron=as.matrix(huron); n=length(huron)

year=seq(1875,by=1,length=n)
par(mfrow=c(2,1))
plot(year,huron,type="o",pch=22,lty=1,pty=2,xaxt="n",xlab="year",ylab="Lake Huron levels")
axis(1,at=seq(1880,1970,by=10),labels=c("1880","1890","1900","1910","1920","1930","1940","1950","1960","1970"))

# Least squares fit
data.huron=data.frame(huron,year); fit.huron=glm(huron~1+year,data=data.huron)
lines(year,fit.huron$fitted.values,lty=1,col="red")

# Residual plot
y=fit.huron$residuals
plot(year,y,type="o",pch=22,lty=1,pty=2,xaxt="n",xlab="year",ylab="",main="Lake Huron residuals from a linear fit")
axis(1,at=seq(1880,1970,by=10),labels=c("1880","1890","1900","1910","1920","1930","1940","1950","1960","1970"))
abline(h=0)

# Sample ACF plot of residuals
acf(y,lag.max=40,ylim=c(-0.6,1),main="Lake Huron residuals from a linear fit")

# Scatter plot of residuals (y_{t}, y_{t+1})
n.res=length(y); plot(y[1:(n.res-1)],y[2:n.res],xlab="y_t",ylab="y_{t+1}")

# Regression of y_{t+1} on y_{t}
res.data.huron=data.frame(y[2:n.res],y[1:(n.res-1)])
res.fit.huron=glm(y[2:n.res]~-1+y[1:(n.res-1)],data=res.data.huron)
(phi.hat=res.fit.huron$coefficients)

# Sample ACF of residuals and confidence bounds
acf.sample=acf(y,lag.max=40,ylim=c(-0.6,1),col="blue",main="Lake Huron residuals from a linear fit")
i=1:40; w.ii=(1-phi.hat^(2*i))*(1+phi.hat^2)/(1-phi.hat^2)-2*i*phi.hat^(2*i)
rho.hat=acf.sample$acf[2:41]
bound.u=rho.hat+1.96*sqrt(w.ii/n); bound.l=rho.hat-1.96*sqrt(w.ii/n)
lines(i,bound.u,type="l",lty=2,pty=2,col="red")
lines(i,bound.l,type="l",lty=2,pty=2,col="red")
lines(i,phi.hat^i,type="p",pch=16) # Model ACF of residuals
legend(20,1,c("sample ACF","95% Conf Bds","model ACF"),col=c("blue","red","black"),lty=c(1,2,0),pch=c(NA_integer_,NA_integer_,16))


