#* Auxiliary functions

source("/sim_algorithm.R")
set.seed(13)

#* Locations on the paraboloid

n1=500;
n2=500;
ang=seq(0,2*pi,l=n1);
t=seq(0,1,l=n2)

coord=mesh(ang,t);
u=coord$x;
v=coord$y;

x=c(v)*cos(c(u));
y=c(v)*sin(c(u));
z=c(v)^2;
sites=cbind(x,y,z);

#* Simulation

df=sim_paraboloid_3d(N=500,tau=c(10,100),nu=c(2,2),s=c(1,1,0.6),sites)

#* Visualization

coordx=v*cos(u);
coordy=v*sin(u);
coordz=v^2;
realization1=array(dim=c(n1,n2),df$component1);
realization2=array(dim=c(n1,n2),df$component2);

surf3D(coordx,coordy,coordz,colvar=realization1,theta=40,phi=40,colkey=T);
plotrgl(lighting=T,smooth=F);
view3d(theta=0,phi=-60,zoom=0.55)

surf3D(coordx,coordy,coordz,colvar=realization2,theta=150,phi=70,colkey=T);
plotrgl(lighting=T,smooth=F);
view3d(theta=0,phi=-60,zoom=0.55)

