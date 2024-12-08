source("/sim_algorithm.R")
set.seed(1234)

ang=seq(0,2*pi,l=100);
t=seq(0,1,l=100)

coord=mesh(ang,t);
u=coord$x;
v=coord$y;

x=c(v)*cos(c(u));
y=c(v)*sin(c(u));
z=c(v)^2;
sites=cbind(x,y,z);

ylab="cumulative distribution function"
xl=c(-3,3)
yl=c(0,1)

nrep=100;
for(j in 1:nrep){
  df=sim_paraboloid_3d(N=100,tau=c(300,300),nu=c(0.5,0.5),s=c(1,1,-0.8),sites)
  e=ecdf(df$component1);
  plot(e,col="lightblue",cex=0.1,xlim=xl,ylim=yl,,main="",ylab=ylab,cex.lab=1.3,cex.axis=1.3);
  par(new=T)
}
curve(pnorm(x),lwd=3,xlim=xl,ylim=yl,main="",ylab=ylab,cex.lab=1.3,cex.axis=1.3);

