rm(list=ls())
graphics.off()

source("utilitary_functions.r")

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

colo=c("red","blue","grey")

#pdf("dominance_diatom_nano_compare_advection_10sp_with_theory.pdf",width=13)
pdf("dominance_diatom_nano_compare_advection_3sp_with_theory.pdf",width=13)
par(mfrow=c(1,2))

#Diatom
#sim_diatom=list(10,11) #Advection, no advection
sim_diatom=list(0,1) #Advection, no advection
lim_max_diatom=25*10^(-4)*10*2

#Nano
#sim_nano=list(12,13)
sim_nano=list(2,3)
lim_max_nano=1.5*10^(-4)*10*2

tot_sim=list(sim_diatom,sim_nano)
tot_max=list(lim_max_diatom,lim_max_nano)

keep_results=matrix(NA,nrow=10,ncol=9)
colnames(keep_results)=c("Abundances","Dlower_adv_micro","Dlower_noadv_micro","Rthreshold_adv_micro","Rthreshold_noadv_micro","Dlower_adv_nano","Dlower_noadv_nano","Rthreshold_adv_nano","Rthreshold_noadv_nano")
corresponding_names=c(expression(r[95~"%"]),expression(D[threshold]))

for(orga in 1:length(tot_sim)){ #Organism: diatom or nano
	for(adv in 1:length(tot_sim[orga][[1]])){ #Advection or no advection
		nb_simu=tot_sim[orga][[1]][[adv]]
	
        	f_tot=read.table(paste("../simulation/lambda_K_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
		colnames(f_tot)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
        	unique_sp=unique(f_tot$sp1)
		if(length(unique_sp)>length(colo)){colo=rainbow(length(unique_sp))}
        	
		f_count=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
		colnames(f_count)=c("species","abundance")


		if(orga==1){
			yl="dominance"
			ml="Micro"
		}else{
			yl=""
			ml="Nano"
		}
		


		if(adv==1){
			plot(1,1,t="n",xlim=range(f_tot$r[f_tot$r>0]),ylim=range(f_tot$dominance,na.rm=T),xlab="r",ylab=yl,log="x",cex.lab=1.5,cex.axis=1.5,axes=F,main=ml)
			axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
			axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
			axis(2,tcl=0.5,cex.axis=1.5) # left
			box()
		}

		if(orga==1){
			leg=c()
			for(l in 1:length(unique_sp)){
				leg=c(leg,paste("Sp=",unique_sp[l],", ",format(100*f_count$abundance[f_count$species==unique_sp[l]]/sum(f_count$abundance),digits=2),"%",sep=""))
			}	
			legend("bottomleft",leg,col=colo[unique_sp+1],pch=1,bty="n",lwd=2,lty=NA)
			mtext("a",side=3,line=1.5,font=2,at=8*10^(-5))
		}else{
			legend("topright",c(expression(U~tau~"/"~3~"="~0),expression(U~tau~"/"~3~"="~0.5)),col="black",pch=c(NA,1),lty=c(1,NA),bty="n",lwd=2)
			mtext("b",side=3,line=1.5,font=2,at=8*10^(-5))
		}

		for (s1 in 1:length(unique_sp)){
        		f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])
		
			if(adv==2){ #No advection
				if(orga==1){
		       			lines(f_plot$r,f_plot$dominance,col=colo[s1])
				}else{
		       			lines(f_plot$r,f_plot$dominance,col=colo[s1])
				}
			}else{
				if(orga==1){
			       		points(f_plot$r,f_plot$dominance,col=colo[s1])
				}else{
		       			points(f_plot$r,f_plot$dominance,col=colo[s1])
				}
			}

		##THEORY
	        #Param simu
        	f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=",dec='.')
	        f_param[,2]=as.numeric(as.character(f_param[,2]))
        	a_lambda=f_param[f_param[,1]=="growth_rate",2]
	        a_Delta=f_param[f_param[,1]=="Delta",2]
        	a_tau=f_param[f_param[,1]=="tau",2]
        	a_tmax=f_param[f_param[,1]=="tmax",2]
        	Utot=f_param[f_param[,1]=="Utot",2]
        	volume=f_param[f_param[,1]=="volume",2]
	        a_Conc=f_count$abundance[f_count$species==unique_sp[s1]]/volume

		Conc_tot=sum(f_count$abundance)/volume

		th_poisson=4/3*pi*f_plot$r^3

		if(Utot>0){
			gamma=1231
			type_line="l"
			th_K=BBM_cdf(f_plot$r,gamma,a_Delta,a_Conc,a_lambda,a_tau,a_tmax)
			th_dominance=a_Conc*(th_K)/(Conc_tot*th_poisson+a_Conc*(th_K-th_poisson))
			lines(f_plot$r,th_dominance,col=colo[s1],lty=2,t=type_line)
		}else{
			gamma=0
			th_K=K_steady_state(f_plot$r,a_lambda,a_Delta,a_Conc,a_tau)
			th_dominance=a_Conc*(th_K)/(Conc_tot*th_poisson+a_Conc*(th_K-th_poisson))
			lines(f_plot$r,th_dominance,col=colo[s1],lwd=2,t="l",lty=3)
			
			th_K=BBM_cdf(f_plot$r,gamma,a_Delta,a_Conc,a_lambda,a_tau,a_tmax)
			th_dominance=a_Conc*(th_K)/(Conc_tot*th_poisson+a_Conc*(th_K-th_poisson))
			points(f_plot$r,th_dominance,col=colo[s1],lty=2,pch="*")
		}


		}

	abline(v=c(tot_max[[orga]][[1]]),lty=2,col=c("black"),lwd=3)

	} #end loop on advection
} #end loop on organism

dev.off()
