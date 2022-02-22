rm(list=ls())
graphics.off()

source("utilitary_functions.r")

colo=c("red","blue","grey")

pdf("dominance_diatom_nano_compare_advection.pdf",width=13)
par(mfrow=c(1,2))

#Diatom
sim_diatom=list(0,1) #Advection, no advection
lim_max_diatom=25*10^(-4)*10*2

#Nano
sim_nano=list(2,3)
lim_max_nano=1.5*10^(-4)*10*2

tot_sim=list(sim_diatom,sim_nano)
tot_max=list(lim_max_diatom,lim_max_nano)

for(orga in 1:length(tot_sim)){ #Organism: diatom or nano
	for(adv in 1:length(tot_sim[orga][[1]])){ #Advection or no advection
		nb_simu=tot_sim[orga][[1]][[adv]]
	
        	f_tot=read.table(paste("../simulation/lambda_K_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
		colnames(f_tot)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
        	unique_sp=unique(f_tot$sp1)

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
			legend("bottomleft",paste("Sp=",unique_sp),col=colo[unique_sp+1],pch=1,bty="n",lwd=2,lty=NA)
			mtext("a",side=3,line=1.5,font=2,at=8*10^(-5))
		}else{
			legend("topright",c(expression(U~tau~"/"~2~"="~0),expression(U~tau~"/"~2~"="~0.5)),col="black",pch=c(NA,1),lty=c(1,NA),bty="n",lwd=2)
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
		}

	abline(v=c(tot_max[[orga]][[1]]),lty=2,col=c("black"),lwd=3)

	} #end loop on advection
} #end loop on organism

dev.off()
