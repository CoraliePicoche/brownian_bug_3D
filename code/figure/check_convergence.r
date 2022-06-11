rm(list=ls())
graphics.off()

source("theoretical_functions.r")

g_steady_state=function(r,lambda,Delta,Conc,tau){
        D=(Delta^2)/(2*tau)
	res=lambda/(4*pi*Conc*D*r)+1
	return(res)
}

K_steady_state=function(r,lambda,Delta,Conc,tau){
        D=(Delta^2)/(2*tau)
	res=lambda*r^2/(2*Conc*D)+(4/3)*pi*r^3
	return(res)
}

pow_r=seq(-4,0,length.out=1000)
r=10^pow_r

tau=0.00028
Delta=6.44437e-05
lambda=1
Conc=10000/1000

seq_tmax_pow=seq(2,10,length.out=50)
seq_tmax=10^seq_tmax_pow
palette_grey=rev(gray.colors(50, start = 0, end = 0.9))

pdf("convergence_wo_advection.pdf",width=10)
par(mfrow=c(1,2))
plot(0.1,0.1,xlim=range(r),ylim=c(0.1,10^6.5),t="n",log="xy",xlab="r",ylab="g(r)")

for(it in 1:length(seq_tmax)){
	t=seq_tmax[it]
	tmp_g=G_theoretical(0,r,lambda,Delta,Conc,tau,Tmax=t)
	lines(r,tmp_g,col=palette_grey[it])
}

tmp_g=G_theoretical(0,r,lambda,Delta,Conc,tau,Tmax=10^3)
lines(r,tmp_g,col="black",lwd=2,lty=2)

tmp_g=G_theoretical(0,r,lambda,Delta,Conc,tau,10^6)
lines(r,tmp_g,col="blue",lwd=2,lty=1)

tmp_g=g_steady_state(r,lambda,Delta,Conc,tau)
lines(r,tmp_g,col="red",lwd=2)


plot(0.1,0.1,xlim=range(r),ylim=c(10^(-4.5),10^4),t="n",log="xy",xlab="r",ylab="K(r)")

for(it in 1:length(seq_tmax)){
	t=seq_tmax[it]
	tmp_K=BBM_cdf(r,0,Delta,Conc,lambda,tau,t)
	lines(r,tmp_K,col=palette_grey[it])
}

tmp_K=BBM_cdf(r,0,Delta,Conc,lambda,tau,10^3)
lines(r,tmp_K,col="black",lwd=2,lty=2)

tmp_K=BBM_cdf(r,0,Delta,Conc,lambda,tau,10^6)
lines(r,tmp_K,col="blue",lwd=2,lty=1)

tmp_K=K_steady_state(r,lambda,Delta,Conc,tau)
lines(r,tmp_K,col="red",lwd=2)

legend("topleft",c("t=10^2","t=10^10","t=10^6","t=1000 (simulations)","steady state"),col=c(palette_grey[1],palette_grey[length(palette_grey)],"blue","black","red"),lwd=2,lty=c(1,1,1,2,1),bty="n")

dev.off()
