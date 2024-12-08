#* Locations on the paraboloid

set.seed(123)
n1=30;
n2=30;
ang=seq(0,2*pi,l=n1);
t=seq(0,1,l=n2)

coord=mesh(ang,t);
u=coord$x;
v=coord$y;

x=c(v)*cos(c(u));
y=c(v)*sin(c(u));
z=c(v)^2;
sites=cbind(x,y,z);

psites=cbind(sites[,1],sites[,2],sqrt(1-sites[,3]));
k=nrow(psites);

#* Convergence rates 

N=1000;
nrep=100; 
trunc=c(28,43,68,106,165,258,403);
qq = array(dim=c(nrep,length(trunc))); 

s1=1;
s2=1;
r12=0.8;
tau1=100;
tau2=100;
nu1=1;
nu2=1;

for(rep in 1:nrep){
	
   matrix1 = array(dim=c(k,N+1));    
   matrix2 = array(dim=c(k,N+1));  
 
   for(n in 0:N){
	
        b11=s1*coef(n,tau1,nu1)/sum_coef(tau1,nu1);
        b22=s2*coef(n,tau2,nu2)/sum_coef(tau2,nu2);
        b12=sqrt(s1*s2)*r12*coef(n,min(tau1,tau2),0.5*(nu1+nu2))/sum_coef(min(tau1,tau2),0.5*(nu1+nu2));
        Bn=cbind(c(b11,b12),c(b12,b22));
        aux=sqrt(2*n+1)*sqrtm(Bn)%*%rnorm(2);
        
	w=rnorm(3);  
        w=w/sqrt(sum(w^2)); 
        dot_prod=psites%*%w;

        matrix1[,(n+1)]=aux[1]*legendre_Pl(n,dot_prod); 
        matrix2[,(n+1)]=aux[2]*legendre_Pl(n,dot_prod); 
                  
   }         
    
   error = c();
   for(j in 1:length(trunc)){

       aux1 = matrix1[,(trunc[j]+1):(N+1)]
       aux2 = matrix2[,(trunc[j]+1):(N+1)]
       error[j] = max((rowSums(aux1))^2+(rowSums(aux2))^2);
	   
   }

   qq[rep,]=error;
}

log_qq=log(sqrt(colMeans(qq)))
lm(log_qq~log(trunc))$coefficients[2]
plot(log(trunc),log_qq,xlim=c(3.32,6),ylim=c(-8.5,0.14),type="b")


