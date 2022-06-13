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
th_poisson=4/3*pi*r^3

seq_tmax_pow=seq(2,8,length.out=50)
seq_tmax=10^seq_tmax_pow
palette_grey=rev(gray.colors(50, start = 0, end = 0.9))


#Diatom
sim=list(0,3) #micro, nano, both without advection

pdf("convergence_wo_advection.pdf",height=7.5,width=5.3)
par(mfcol=c(3,2))

for (i in 1:length(sim)){

	if(i==1){
		yl="g(r)"
		tl="Micro"
		let=c("a","c","e")
	}else{
		yl=""
		tl="Nano"
		let=c("b","d","f")
	}

	nb_simu=sim[i]

	#Param simu
	f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=",dec='.')
	f_param[,2]=as.numeric(as.character(f_param[,2]))
	a_lambda=f_param[f_param[,1]=="growth_rate",2]
	a_Delta=f_param[f_param[,1]=="Delta",2]
	a_tau=f_param[f_param[,1]=="tau",2]
	a_lim=f_param[f_param[,1]=="radius",2]*10^2 #Radius is converted to cm here
	lim=a_lim*10*2
	volume=f_param[f_param[,1]=="volume",2]
	a_Conc=f_param[f_param[,1]=="init_size 0",2]/volume


	par(mar=c(2,4,4,2))

	tmp_g=g_steady_state(r,a_lambda,a_Delta,a_Conc,a_tau)
	plot(0.1,0.1,xlim=range(r),ylim=range(tmp_g),t="n",log="xy",xlab="",ylab=yl,main=tl)
	mtext(let[1],side=3,line=0.5,font=2,at=8*10^(-5),cex=0.75)

	lines(r,tmp_g,col="red",lwd=2)

for(it in 1:length(seq_tmax)){
	t=seq_tmax[it]
	tmp_g=G_theoretical(0,r,a_lambda,a_Delta,a_Conc,a_tau,Tmax=t)
	lines(r,tmp_g,col=palette_grey[it])
}

tmp_g=G_theoretical(0,r,a_lambda,a_Delta,a_Conc,a_tau,Tmax=10^3)
lines(r,tmp_g,col="black",lwd=2,lty=2)

tmp_g=G_theoretical(0,r,a_lambda,a_Delta,a_Conc,a_tau,25000) #25000 IS A WEEK
lines(r,tmp_g,col="blue",lwd=2,lty=1)

abline(v=lim,col="orange",lty=3,lwd=3)


	if(i==1){
		yl="K(r)"
	}else{
		yl=""
	}


par(mar=c(3,4,3,2))

tmp_K=K_steady_state(r,a_lambda,a_Delta,a_Conc,a_tau)
plot(0.1,0.1,xlim=range(r),ylim=range(tmp_K),t="n",log="xy",xlab="",ylab=yl)
mtext(let[2],side=3,line=0.5,font=2,at=8*10^(-5),cex=0.75)

lines(r,tmp_K,col="red",lwd=2)

for(it in 1:length(seq_tmax)){
        t=seq_tmax[it]
        tmp_K=BBM_cdf(r,0,a_Delta,a_Conc,a_lambda,a_tau,t)
        lines(r,tmp_K,col=palette_grey[it])
}

tmp_K=BBM_cdf(r,0,a_Delta,a_Conc,a_lambda,a_tau,10^3)
lines(r,tmp_K,col="black",lwd=2,lty=2)

tmp_K=BBM_cdf(r,0,a_Delta,a_Conc,a_lambda,a_tau,25000)
lines(r,tmp_K,col="blue",lwd=2,lty=1)

abline(v=lim,col="orange",lty=3,lwd=3)
	

	if(i==1){
		yl="dominance(r)"
	}else{
		yl=""
	}


par(mar=c(4,4,2,2))
plot(0.1,0.1,xlim=range(r),ylim=c(0.3,1),t="n",log="x",xlab="r",ylab=yl)
mtext(let[3],side=3,line=0.5,font=2,at=8*10^(-5),cex=0.75)

for(it in 1:length(seq_tmax)){
	t=seq_tmax[it]
	th_K=BBM_cdf(r,0,a_Delta,a_Conc,a_lambda,a_tau,t)
	th_dominance=a_Conc*(th_K)/(a_Conc*3*th_poisson+a_Conc*(th_K-th_poisson))
	lines(r,th_dominance,col=palette_grey[it])
}

th_K=BBM_cdf(r,0,a_Delta,a_Conc,a_lambda,a_tau,10^3)
th_dominance=a_Conc*(th_K)/(a_Conc*3*th_poisson+a_Conc*(th_K-th_poisson))
lines(r,th_dominance,col="black",lwd=2,lty=2)

th_K=BBM_cdf(r,0,a_Delta,a_Conc,a_lambda,a_tau,25000)
th_dominance=a_Conc*(th_K)/(a_Conc*3*th_poisson+a_Conc*(th_K-th_poisson))
lines(r,th_dominance,col="blue",lwd=2,lty=1)

th_K=K_steady_state(r,a_lambda,a_Delta,a_Conc,a_tau)
th_dominance=a_Conc*(th_K)/(a_Conc*3*th_poisson+a_Conc*(th_K-th_poisson))
lines(r,th_dominance,col="red",lwd=2)

abline(v=lim,col="orange",lty=3,lwd=3)

legend("bottomleft",c(expression(t=10^2),expression(t=10^9),"t=25000","t=1000","steady state"),col=c(palette_grey[1],palette_grey[length(palette_grey)],"blue","black","red"),lwd=2,lty=c(1,1,1,2,1),bty="n")

}

dev.off()

