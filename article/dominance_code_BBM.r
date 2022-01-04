rm(list=ls())
graphics.off()

library(pracma) #erf function

BBM_cdf=function(r,gamma,D,Conc,lambda,t=NA){
	res=NA
	A=4/3*pi*r^3
	if(gamma>0){
		B=r^2/(6*D)
		C=(sqrt(gamma)/(6*sqrt(2)*(D^1.5)))*(r^3)*atan(r*sqrt(gamma/(2*D)))
		Dd=log((gamma*r^2)/(2*D)+1)/(6*gamma)
		E=(sqrt(gamma)*pi*r^3)/(12*sqrt(2)*D^1.5)
		res=A+(2*lambda/Conc)*(B+C+Dd-E)
	}else if(gamma==0){
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

#Parameters Diatoms
a_lambda=1*0.0002
a_C=10000/1000
a_r=25*10^(-6)

#Parameter Nano
#a_lambda=2.5*0.0002
#a_C=10000/10
#a_r=1.5*10^(-6)

#Common parameters
a_D=8.314*293/(6.0225*10^23)*1/(6*pi*a_r*10^(-3))
n_species=3
seq_r=seq(10^(-4),10^(-3),length.out=100)
seq_r=c(seq_r,seq(10^(-3),10^(-2),length.out=100))
seq_r=c(seq_r,seq(10^(-2),10^(-1),length.out=100))
seq_r=c(seq_r,seq(10^(-1),1,length.out=100))

concentration=rep(a_C,n_species)

dom=array(NA,dim=c(n_species,length(seq_r),2))

den_1=sum(concentration)*poisson_cdf(seq_r)

for (i in 1:n_species){
	num=concentration[i]*(BBM_cdf(seq_r,0.5,a_D,concentration[i],a_lambda))
	dom[i,,1]=num/(den_1+num-poisson_cdf(seq_r)*concentration[i])
	num=concentration[i]*(BBM_cdf(seq_r,0.0,a_D,concentration[i],a_lambda,1000))
	dom[i,,2]=num/(den_1+num-poisson_cdf(seq_r)*concentration[i])
}

pdf("dominance_BBM_diat.pdf",width=10)
par(mfrow=c(2,2))
plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n",log="x",main="Advection")
points(seq_r,dom[1,,1],lty=1,col="black",pch=16)
points(seq_r,dom[2,,1],lty=2,col="blue",pch=16)
points(seq_r,dom[3,,1],lty=3,col="red",pch=16)
abline(h=1/3,lty=2,col="grey",lwd=2)
legend("bottomleft",c("sp1=33%","sp2=33%","sp3=33%"),pch=16,col=c("black","blue","red"),bty="n")

plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n",main="No advection",log="x")
points(seq_r,dom[1,,2],lty=1,col="black",pch=16)
points(seq_r,dom[2,,2],lty=2,col="blue",pch=16)
points(seq_r,dom[3,,2],lty=3,col="red",pch=16)
abline(h=1/3,lty=2,col="grey",lwd=2)

mtext("Nano", outer = TRUE, cex = 1.5,line=-1.75)

#"Skewed distributions"
concentration=c(1000,3333,5667)/1000 

dom=array(NA,dim=c(n_species,length(seq_r),2))

den_1=sum(concentration)*poisson_cdf(seq_r)

for (i in 1:n_species){
        num=concentration[i]*(BBM_cdf(seq_r,0.5,a_D,concentration[i],a_lambda))
        dom[i,,1]=num/(den_1+num-poisson_cdf(seq_r)*concentration[i])
        num=concentration[i]*(BBM_cdf(seq_r,0.0,a_D,concentration[i],a_lambda,1000))
        dom[i,,2]=num/(den_1+num-poisson_cdf(seq_r)*concentration[i])
}

plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n",log="x",main="Advection")
points(seq_r,dom[1,,1],lty=1,col="black",pch=16)
points(seq_r,dom[2,,1],lty=2,col="blue",pch=16)
points(seq_r,dom[3,,1],lty=3,col="red",pch=16)
legend("bottomleft",c("sp1=10%","sp2=33%","sp3=57%"),pch=16,col=c("black","blue","red"),bty="n")

plot(0,0,xlim=range(seq_r),ylim=c(0,1),xlab="r",ylab="dominance",t="n",main="No advection",log="x")
points(seq_r,dom[1,,2],lty=1,col="black",pch=16)
points(seq_r,dom[2,,2],lty=2,col="blue",pch=16)
points(seq_r,dom[3,,2],lty=3,col="red",pch=16)


dev.off()
