rm(list=ls())
graphics.off()

source("utilitary_functions.r")

colo=c("red","blue","grey")

pdf("dominance_diatom_nano_compare_advection_10sp.pdf",width=13)
#pdf("dominance_diatom_nano_compare_advection.pdf",width=13)
par(mfrow=c(1,2))

#Diatom
sim_diatom=list(10,11) #Advection, no advection
#sim_diatom=list(0,1) #Advection, no advection
lim_max_diatom=25*10^(-4)*10*2

#Nano
sim_nano=list(12,13)
#sim_nano=list(2,3)
lim_max_nano=1.5*10^(-4)*10*2

tot_sim=list(sim_diatom,sim_nano)
tot_max=list(lim_max_diatom,lim_max_nano)

keep_results=matrix(NA,nrow=10,ncol=9)
colnames(keep_results)=c("Abundances","Dlower_adv_micro","Dlower_noadv_micro","Rthreshold_adv_micro","Rthreshold_noadv_micro","Dlower_adv_nano","Dlower_noadv_nano","Rthreshold_adv_nano","Rthreshold_noadv_nano")
corresponding_names=c("r/D_i<95%","D_i/r=10 diam")

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
		
			#SOME ANALYSES
			id_1=which(f_plot$dominance<0.95)[1]
			keep_results[s1,1+adv+4*(orga-1)]=f_plot$r[id_1]
			id_2=1
			u_1=f_plot$r[id_2]
			while (u_1<tot_max[[orga]][[1]]){
				id_2=id_2+1
				u_1=f_plot$r[id_2]
			}
			keep_results[s1,1+2+adv+4*(orga-1)]=f_plot$dominance[id_2]
			
			keep_results[s1,1]=f_count$abundance[f_count$species==unique_sp[s1]]
		
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
write.table(keep_results,"dom_indices_10sp.txt",sep=";",dec=".",col.names=TRUE,row.names=F)

###Part where we use kee_results
pdf("carac_10species.r")
par(mfrow=c(2,2),mar=c(2,4,4,2))
plot(0,0,t="n",xlim=range(keep_results[,"Abundances"]),ylim=range(keep_results[,c("Dlower_adv_micro","Dlower_noadv_micro")]),xlab="",ylab=corresponding_names[1],main="Micro",log="x")
points(keep_results[,"Abundances"],keep_results[,"Dlower_adv_micro"],pch=16)
lines(keep_results[,"Abundances"],keep_results[,"Dlower_noadv_micro"],lty=1)
mtext("a",side=3,line=1.5,font=2,at=8*10^(-5))

plot(0,0,t="n",xlim=range(keep_results[,"Abundances"]),ylim=range(keep_results[,c("Dlower_adv_nano","Dlower_noadv_nano")]),xlab="",ylab="",main="Nano",log="x")
points(keep_results[,"Abundances"],keep_results[,"Dlower_adv_nano"],pch=16)
lines(keep_results[,"Abundances"],keep_results[,"Dlower_noadv_nano"],lty=1)
mtext("b",side=3,line=1.5,font=2,at=8*10^(-5))

par(mar=c(4,4,2,2))
plot(0,0,t="n",xlim=range(keep_results[,"Abundances"]),ylim=range(keep_results[,c("Rthreshold_adv_micro","Rthreshold_noadv_micro","Rthreshold_adv_nano","Rthreshold_noadv_nano")]),xlab="abundances",ylab=corresponding_names[2],log="x")
points(keep_results[,"Abundances"],keep_results[,"Rthreshold_adv_micro"],pch=16)
lines(keep_results[,"Abundances"],keep_results[,"Rthreshold_noadv_micro"],lty=1)
mtext("c",side=3,line=0.5,font=2,at=8*10^(-5))

plot(0,0,t="n",xlim=range(keep_results[,"Abundances"]),ylim=range(keep_results[,c("Rthreshold_adv_micro","Rthreshold_noadv_micro","Rthreshold_adv_nano","Rthreshold_noadv_nano")]),xlab="abundances",ylab="",log="x")
points(keep_results[,"Abundances"],keep_results[,"Rthreshold_adv_nano"],pch=16)
lines(keep_results[,"Abundances"],keep_results[,"Rthreshold_noadv_nano"],lty=1)
mtext("d",side=3,line=0.5,font=2,at=8*10^(-5))
legend("bottomright",c("Adv","No adv"),pch=c(16,NA),lty=c(NA,1))
dev.off()
