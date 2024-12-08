source("/sim_algorithm.R")
set.seed(13)

time=c();
n=c(100,200,300,400,500);
N=100;

for(j in 1:length(n)){
  
   ang=seq(0,2*pi,l=n[j]);
   t=seq(0,1,l=n[j])

   coord=mesh(ang,t);
   u=coord$x;
   v=coord$y;

   x=c(v)*cos(c(u));
   y=c(v)*sin(c(u));
   z=c(v)^2;
   sites=cbind(x,y,z);

   t1=Sys.time()
     df=sim_paraboloid_3d(N=N,tau=c(10,100),nu=c(2,2),s=c(1,1,0.6),sites)
   t2=Sys.time()

   time[j]=t2-t1; 
   print(j)
}
