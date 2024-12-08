#* Libraries

require("plot3D");
require("gsl");
require("expm");
require("plot3Drgl")

#* Coefficients circular-Matérn model

coef<-function(n,tau,nu){(n^2+tau)^(-nu-0.5)}
sum_coef<-function(tau,nu){sum(((0:1000)^2+tau)^(-nu-0.5))}

#* Simulation algorithm paraboloid

sim_paraboloid_3d <- function(N,tau,nu,s,sites){
   
   psites=cbind(sites[,1],sites[,2],sqrt(1-sites[,3]));
   component1=rep(0,nrow(psites));
   component2=rep(0,nrow(psites));
   
   for(n in 0:N){
	
     b11=s[1]*coef(n,tau[1],nu[1])/sum_coef(tau[1],nu[1]);
     b22=s[2]*coef(n,tau[2],nu[2])/sum_coef(tau[2],nu[2]);
     b12=sqrt(s[1]*s[2])*s[3]*coef(n,min(tau[1],tau[2]),0.5*(nu[1]+nu[2])) / sum_coef(min(tau[1],tau[2]),0.5*(nu[1]+nu[2]));
     Bn=cbind(c(b11,b12),c(b12,b22));
     aux=sqrt(2*n+1)*sqrtm(Bn)%*%rnorm(2);
     w=rnorm(3);  
     w=w/sqrt(sum(w^2)); 
     dot_prod=psites%*%w;
  
     component1=component1+aux[1]*legendre_Pl(n,dot_prod); 
     component2=component2+aux[2]*legendre_Pl(n,dot_prod); 
                  
   }
   
   return(data.frame(sites,component1,component2))
}

#* Simulation algorithm parabola

sim_parabola_2d <- function(N,tau,nu,s,sites){
   
   psites=cbind(sites[,1],sqrt(1-sites[,2]));
   component1=rep(0,nrow(psites));
   component2=rep(0,nrow(psites));
   
   for(n in 0:N){
	
     b11=s[1]*coef(n,tau[1],nu[1])/sum_coef(tau[1],nu[1]);
     b22=s[2]*coef(n,tau[2],nu[2])/sum_coef(tau[2],nu[2]);
     b12=sqrt(s[1]*s[2])*s[3]*coef(n,min(tau[1],tau[2]),0.5*(nu[1]+nu[2])) / sum_coef(min(tau[1],tau[2]),0.5*(nu[1]+nu[2]));

     Bn=cbind(c(b11,b12),c(b12,b22));
     aux=ifelse(n>0,sqrt(2),1)*sqrtm(Bn)%*%rnorm(2);
     w=runif(1,0,2*pi);  
     dot_prod=psites%*%c(cos(w),sin(w));
  
     component1=component1+aux[1]*cos(n*acos(dot_prod)); 
     component2=component2+aux[2]*cos(n*acos(dot_prod)); 
                  
   }
   
   return(data.frame(sites,component1,component2))
}
