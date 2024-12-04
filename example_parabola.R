#* Auxiliary functions

source("/Users/alfredo/Desktop/paraboloid/new/final/sim_algorithm.R")
library("latticeExtra")
library("sp")
set.seed(123)

#* Locations on the parabola

n1=1000;
t=seq(-1,1,l=n1)

x=t;
y=x^2;
sites=cbind(x,y);

#* Simulation

df=sim_parabola_2d(N=500,tau=c(100,100),nu=c(0.5,1.5),s=c(1,1,0.65),sites)
coordinates(df) <- ~x+y

#* Visualization

p<-spplot(df,colorkey=T,xlab="coordx",ylab="coordy",
             ylim=c(-0.6,1.6),cex=8,pch="|")

p+layer(panel.points(x,y,type="l",lwd=3),data=df)

