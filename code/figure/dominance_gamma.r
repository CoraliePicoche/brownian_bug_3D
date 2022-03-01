rm(list=ls())
graphics.off()

source("utilitary_functions.r")
source("theoretical_functions.r")

pdf("theoretical_dominance_with_adv.pdf",width=13)
par(mfrow=c(1,2))

#Diatom
sim_diatom=list(0,1) #Advection, no advection
lim_max_diatom=25*10^(-4)*10*2

#Nano
sim_nano=list(2,3)
lim_max_nano=1.5*10^(-4)*10*2

tot_sim=list(sim_diatom,sim_nano)
tot_max=list(lim_max_diatom,lim_max_nano)

u_r=seq(10^(-4),1,length.out=1000)
th_poisson=4/3*pi*u_r^3

for(orga in 1:length(tot_sim)){ #Organism: diatom or nano
        for(adv in 1:length(tot_sim[orga][[1]])){ #Advection or no advection
		nb_simu=tot_sim[orga][[1]][[adv]]
		f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),sep="=",header=F,dec=".")
		f_param[,2]=as.numeric(as.character(f_param[,2]))
		colnames(f_param)=c("name","value")
		volume=f_param[f_param$name=="volume","value"]

		Utot=f_param[f_param$name=="Utot","value"]
if(Utot==0.5){
        gamma=1231
}else if (Utot==0.0){
        gamma=0
}
a_Delta=f_param[f_param$name=="Delta","value"]
a_lambda=f_param[f_param$name=="growth_rate","value"]
a_tau=f_param[f_param$name=="tau","value"]
a_tmax=f_param[f_param$name=="tmax","value"]
                
if(orga==1){
                        yl="dominance"
                        ml="Micro"
			concentration=10
                }else{
                        yl=""
                        ml="Nano"
			concentration=10^3
                }

th_bbm=BBM_cdf(u_r,gamma,a_Delta,concentration,a_lambda,a_tau,t=a_tmax)
th_dominance=concentration*th_bbm/(3*concentration*th_poisson+concentration*(th_bbm-th_poisson))

if(Utot==0.5){
	gamma_bis=164
}else if (Utot==0.0){
        gamma_bis=0
}
tau_bis=0.0021
th_bbm_bis=BBM_cdf(u_r,gamma_bis,a_Delta,concentration,a_lambda,tau_bis,t=a_tmax)
th_dominance_bis=concentration*th_bbm_bis/(3*concentration*th_poisson+concentration*(th_bbm_bis-th_poisson))

                if(adv==1){
                        plot(1,1,t="n",xlim=range(u_r),ylim=c(0.3,1),xlab="r",ylab=yl,log="x",cex.lab=1.5,cex.axis=1.5,axes=F,main=ml)
                        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
                        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
                        axis(2,tcl=0.5,cex.axis=1.5) # left
                        box()
                }

                if(orga==1){
                        legend("bottomleft",c(eval(substitute(expression(paste(tau,"=",a_tau,sep="")),list(a_tau=a_tau))),eval(substitute(expression(paste(tau,"=",a_tau,sep="")),list(a_tau=tau_bis)))),col=c("black","grey"),lty=c(1),bty="n",lwd=2)
                        mtext("a",side=3,line=1.5,font=2,at=8*10^(-5))
                }else{
                        legend("topright",c(expression(U~tau~"/"~3~"="~0),expression(U~tau~"/"~3~"="~0.5)),col="black",pch=c(NA,1),lty=c(1,NA),bty="n",lwd=2)
                        mtext("b",side=3,line=1.5,font=2,at=8*10^(-5))
                }

                        if(adv==2){ #No advection
                             lines(u_r,th_dominance,col="black")
                             lines(u_r,th_dominance_bis,col="grey")
                        }else{
                             points(u_r,th_dominance,col="black")
                             points(u_r,th_dominance_bis,col="grey")
                        }

        abline(v=c(tot_max[[orga]][[1]]),lty=2,col=c("black"),lwd=3)

        } #end loop on advection
} #end loop on organism

dev.off()
