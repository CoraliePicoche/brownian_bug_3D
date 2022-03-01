rm(list=ls())
graphics.off()

source("utilitary_functions.r")
source("theoretical_functions.r")

nb_simu_tot=c(1,3)
palette_grey=rev(gray.colors(100, start = 0, end = 0.75))

s=2 #We choose only one species
pow_th=seq(-4,0,length.out=50)
th_r=10^pow_th
th_poisson=4/3*pi*th_r^3

pdf("theory_dominance.pdf",width=10)
par(mfcol=c(1,2))

for(n in 1:length(nb_simu_tot)){
	    nb_simu=nb_simu_tot[n]

if(n==1){
	yl="dominance"
 	ml="Micro"
}else{
      	yl=""
       	ml="Nano"
}


#Let's add the dominance index for Thomas distribution
f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),sep="=",header=F,dec=".")
colnames(f_param)=c("name","value")
f_param[,2]=as.numeric(as.character(f_param[,2]))
gamma=0
a_Delta=f_param[f_param$name=="Delta","value"]
a_lambda=f_param[f_param$name=="growth_rate","value"]
a_tau=f_param[f_param$name=="tau","value"]
a_volume=f_param[f_param$name=="volume","value"]
a_radius=f_param[f_param$name=="radius","value"]

concentrations=c(f_param[f_param$name=="init_size 0","value"],f_param[f_param$name=="init_size 1","value"],f_param[f_param$name=="init_size 2","value"])/a_volume

pow_time=seq(2,5,length.out=100)
seq_time=10^pow_time

plot(0,0,xlim=range(th_r),ylim=c(0.33,1),ylab=yl,xlab="r",t="n",log="x",main=ml,cex.main=1.25)
for(time in 1:length(seq_time)){
        th_bbm=BBM_cdf(th_r,gamma,a_Delta,concentrations[s],a_lambda,a_tau,t=seq_time[time])
        th_dominance=concentrations[s]*(th_bbm)/(sum(concentrations)*th_poisson+concentrations[s]*(th_bbm-th_poisson))
        lines(th_r,th_dominance,col=palette_grey[time])
}

gamma=1231
th_bbm=BBM_cdf(th_r,gamma,a_Delta,concentrations[s],a_lambda,a_tau)
th_dominance=concentrations[s]*(th_bbm)/(sum(concentrations)*th_poisson+concentrations[s]*(th_bbm-th_poisson))
lines(th_r,th_dominance,col="red")
abline(v=a_radius*2*10*10^2,lty=2,lwd=2)

if(n==1){
 legend("bottomleft",c("No adv t=100","No adv t=100000","Advection"),col=c(palette_grey[1],palette_grey[100],"red"),bty="n",lwd=3,lty=1,pch=NA)
 mtext("a",side=3,line=1.5,font=2,at=8*10^(-5),cex=1.25)
}else{
 mtext("b",side=3,line=1.5,font=2,at=8*10^(-5),cex=1.25)
}

}

dev.off()
