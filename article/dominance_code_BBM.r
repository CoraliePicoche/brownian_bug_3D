rm(list=ls())
graphics.off()

library(pracma) #erf function


G_theoretical=function(r,gamma,lambda,Delta,C_0,tau,Tmax=NA){ #U=0 in the absence of advection
  D=(Delta^2)/(2*tau)
  Tmax=Tmax*tau

  #a_r=25*10^(-6)
  #D=8.314*293/(6.0225*10^23)*1/(6*pi*a_r*10^(-3)) 

  if(gamma>0){
    tmp=lambda/(2*pi*C_0)*(sqrt(gamma)*atan(sqrt(gamma)*r/sqrt(2*D))/(2^1.5*D^1.5)+1/(2*D*r)-pi*sqrt(gamma)/(2^2.5*D^1.5))+1
  }else{ #This corresponds to the case U=0
        if(is.na(Tmax)){
                print("Tmax should have a value")
                stop()
        }
    tmp=lambda/(4*D*pi*C_0*r)*(1-erf(r/(sqrt(8*Tmax*D))))+1
  }
  return(tmp)
}

G_theoretical_4pi=function(r,gamma,lambda,Delta,C_0,tau,Tmax=NA){ #U=0 in the absence of advection
  D=(Delta^2)/(2*tau)
  Tmax=Tmax*tau

  if(gamma>0){
    tmp=lambda/(2*pi*C_0)*(sqrt(gamma)*atan(sqrt(gamma)*r/sqrt(2*D))/(2^1.5*D^1.5)+1/(2*D*r)-pi*sqrt(gamma)/(2^2.5*D^1.5))+1
  }else{ #This corresponds to the case U=0
        if(is.na(Tmax)){
                print("Tmax should have a value")
                stop()
        }
    tmp=lambda/(4*D*pi*C_0*r)*(1-erf(r/(sqrt(8*Tmax*D))))+1
  }
  tmp=4*pi*r^2*tmp
  return(tmp)
}


BBM_cdf=function(r,gamma,Delta,Conc,lambda,tau,t=NA){
        D=(Delta^2)/(2*tau)
	res=NA
	A=4/3*pi*r^3
  	t=t*tau

	#a_r=25*10^(-6)
	#D=8.314*293/(6.0225*10^23)*1/(6*pi*a_r*10^(-3))

	if(gamma>0){
		B=r^2/(6*D)
		C=(sqrt(gamma)/(6*sqrt(2)*(D^1.5)))*(r^3)*atan(r*sqrt(gamma/(2*D)))
		Dd=log((gamma*r^2)/(2*D)+1)/(6*gamma)
		E=(sqrt(gamma)*pi*r^3)/(12*sqrt(2)*D^1.5)
		res=A+(2*lambda/Conc)*(B+C+Dd-E)
	}else{
		B=r^2/2
		C=1/2*erf(r/(sqrt(8*D*t)))*(r^2-4*D*t)
		Dd=sqrt(2*D*t)*r/sqrt(pi)*exp(-r^2/(8*D*t))
		res=A+(lambda/(Conc*D))*(B-C-Dd)
	}
	return(res)
}

poisson_cdf=function(r){
	a=4/3*pi*r^3
	return(a)
}

#Parameters Diatoms / Nano
a_Delta=6.44437e-05
#a_Delta=0.00026309
lim=25*10^(-4)*10
#lim=1.5*10^(-4)*10
a_lambda=1
#a_lambda=2.5

a_tau=0.00028

#Test to compute D
a_gamma=1231

a_C=10000/1000

#Common parameters
n_species=3
seq_r=seq(10^(-4),10^(-3),length.out=100)
seq_r=c(seq_r,seq(10^(-3),10^(-2),length.out=100))
seq_r=c(seq_r,seq(10^(-2),10^(-1),length.out=100))
seq_r=c(seq_r,seq(10^(-1),1,length.out=100))

concentration=rep(a_C,n_species)

dom=array(NA,dim=c(n_species,length(seq_r),2))

den_1=sum(concentration)*poisson_cdf(seq_r)

pdf("K_dominance_BBM_r_micro.pdf",width=12)
par(mfcol=c(2,3),oma=c(1,1,1,1),cex=0.85)
#plot(seq_r,G_theoretical(seq_r,a_gamma,a_lambda,a_Delta,concentration[1],a_tau),log="xy",ylab="pcf",xlab="r",t="l")
plot(seq_r,BBM_cdf(seq_r,a_gamma,a_Delta,concentration[1],a_lambda,a_tau),t="l",xlab="r",ylab="K",log="xy",main="Advection")
plot(seq_r,BBM_cdf(seq_r,0,a_Delta,concentration[1],a_lambda,a_tau,t=1000),t="l",xlab="r",ylab="K",log="xy",main="No Advection")
#points(seq_r,piou,col="red",pch=1,cex=0.5)
#legend("topleft",c("Theory","R integrate"),col=c("black","red"),lty=c(1,NA),pch=c(NA,1),bty="n")

for(i in 1:n_species){
	num=concentration[i]*(BBM_cdf(seq_r,a_gamma,a_Delta,concentration[i],a_lambda,a_tau))
	dom[i,,1]=num/(den_1+num-poisson_cdf(seq_r)*concentration[i])
	num=concentration[i]*(BBM_cdf(seq_r,0.0,a_Delta,concentration[i],a_lambda,a_tau,1000))
	dom[i,,2]=num/(den_1+num-poisson_cdf(seq_r)*concentration[i])
}

plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n",log="x",main="Advection")
points(seq_r,dom[1,,1],lty=1,col="red",pch=16)
points(seq_r,dom[2,,1],lty=2,col="blue",pch=16)
points(seq_r,dom[3,,1],lty=3,col="black",pch=16)
abline(h=1/3,lty=2,col="grey",lwd=2)
legend("bottomleft",c("sp1=33%","sp2=33%","sp3=33%"),pch=16,col=c("red","blue","black"),bty="n")
abline(v=lim,col="grey",lty=3,lwd=3)

plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n",main="No advection",log="x")
points(seq_r,dom[1,,2],lty=1,col="red",pch=16)
points(seq_r,dom[2,,2],lty=2,col="blue",pch=16)
points(seq_r,dom[3,,2],lty=3,col="black",pch=16)
abline(h=1/3,lty=2,col="grey",lwd=2)
abline(v=lim,col="grey",lty=3,lwd=3)

mtext("Micro", outer = TRUE, cex = 1.5,line=-1.5)

#"Skewed distributions"
concentration=c(1000,3333,5667)/1000

dom=array(NA,dim=c(n_species,length(seq_r),2))

den_1=sum(concentration)*poisson_cdf(seq_r)

for (i in 1:n_species){
        num=concentration[i]*(BBM_cdf(seq_r,a_gamma,a_Delta,concentration[i],a_lambda,a_tau))
        dom[i,,1]=num/(den_1+num-poisson_cdf(seq_r)*concentration[i])
        num=concentration[i]*(BBM_cdf(seq_r,0.0,a_Delta,concentration[i],a_lambda,a_tau,1000))
        dom[i,,2]=num/(den_1+num-poisson_cdf(seq_r)*concentration[i])
}

plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n",log="x",main="Advection")
points(seq_r,dom[1,,1],lty=1,col="red",pch=16)
points(seq_r,dom[2,,1],lty=2,col="blue",pch=16)
points(seq_r,dom[3,,1],lty=3,col="black",pch=16)
legend("bottomleft",c("sp1=10%","sp2=33%","sp3=57%"),pch=16,col=c("red","blue","black"),bty="n")
abline(v=lim,col="grey",lty=3,lwd=3)

plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n",main="No advection",log="x")
points(seq_r,dom[1,,2],lty=1,col="red",pch=16)
points(seq_r,dom[2,,2],lty=2,col="blue",pch=16)
points(seq_r,dom[3,,2],lty=3,col="black",pch=16)
abline(v=lim,col="grey",lty=3,lwd=3)

dev.off()
