rm(list=ls())
graphics.off()

library(spatstat)
source("utilitary_functions.r")
source("theoretical_functions.r")

colo=c("red","blue","black")

nb_simu=150

a_gamma=1231

show_dist=FALSE

pdf("compare_BBM_spatstat_bandwidth_whole_range.pdf",width=15,height=5)
par(mfrow=c(1,3),mar=c(4,5,3,4))

f=read.table(paste("../simulation/lambda_K_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
unique_sp=unique(f$sp1)
if (length(unique_sp)>length(colo)){
        colo=rainbow(length(unique_sp))
}

#Count species
f_count=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";",dec='.')
colnames(f_count)=c("species","abundance")

#Param simu
f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=",dec='.')
f_param[,2]=as.numeric(as.character(f_param[,2]))
delta_val=format(f_param[f_param[,1]=="delta",2],digits=2)
a_lambda=f_param[f_param[,1]=="growth_rate",2]
a_Delta=f_param[f_param[,1]=="Delta",2]
a_tau=f_param[f_param[,1]=="tau",2]

#Spatial distribution
tab_repart=read.table(paste("../simulation/Spatial_distribution_",nb_simu,".txt",sep=""), header=F, sep=";")
colnames(tab_repart)=c("t","x","y","z","yfirst","ancestor","species")

Concentration=f_count$abundance/f_param[f_param[,1]=="volume",2]

max_r=max(f$r)

for (sp in unique(f_count$species)){
	xl="r"
	if(sp==0){
		yl="PCF"
	}else{
		yl=""
	}
	f_plot=subset(f,sp1==sp&sp2==sp)
  	
	tab_repart_tmp=subset(tab_repart,species==sp)
    	ppp_repart=pp3(tab_repart_tmp$x,y=tab_repart_tmp$y,z=tab_repart_tmp$z,box3(c(0,1)))
	
	if(show_dist){
		dist=unique(c(pairdist(ppp_repart)))
		a=sort(dist[dist>0],decreasing=F)
		a=a[1:1000]
		seq_r=unique(f_plot$r[f_plot$r<10^(-2)])
		tab_dist=rep(0,length(seq_r))
		for (i in 1:length(a)){
			if (a[i]<=seq_r[1]){
				tab_dist[1]=tab_dist[1]+1
	       		}else{
				j=2
				found=FALSE
				while (j<=length(seq_r) && !found){
					found=(a[i]>seq_r[j-1] && a[i] <=seq_r[j])
					if (found) tab_dist[j]=tab_dist[j]+1
					j=j+1
       				}
			}
		}
	}
#PCF
  	pcf_spat_3=pcf3est(ppp_repart,rmax=max_r,nrval=5000,delta=10^(-3))
  	pcf_spat_3p5=pcf3est(ppp_repart,rmax=max_r,nrval=5000,delta=10^(-3.2))
  	pcf_spat_4=pcf3est(ppp_repart,rmax=max_r,nrval=5000,delta=10^(-4))
  	pcf_spat_5=pcf3est(ppp_repart,rmax=max_r,nrval=5000,delta=10^(-5))
 # 	pcf_spat_6=pcf3est(ppp_repart,rmax=10^(-2),nrval=1000)

	#Theory
	th_bbm=G_theoretical(a_gamma,f_plot$r,a_lambda,a_Delta,Concentration[sp+1],a_tau)

	plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab=yl,axes=F,xlim=c(10^(-4),max_r),cex.lab=1.5,ylim=c(range(range(th_bbm),range(f$pcf[f$pcf>0]))),main=paste("Species =",sp))
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()
#	points(pcf_spat_3$r,pcf_spat_3$trans,col="orchid",cex=1.5,pch=16)
	points(pcf_spat_3p5$r,pcf_spat_3p5$trans,col="orchid",cex=1.5,pch=16)
	points(pcf_spat_4$r,pcf_spat_4$trans,col="purple",cex=1.5,pch=17)
	points(pcf_spat_5$r,pcf_spat_5$trans,col="orange",cex=1.5,pch=18)
#	points(pcf_spat_6$r,pcf_spat_6$trans,col="chartreuse",cex=1.5,pch=19)
        lines(f_plot$r,f_plot$pcf,col=colo[sp+1],lwd=2)
	lines(f_plot$r,th_bbm,lty=2,col="grey",lwd=2)

	if(show_dist){
		par(new = TRUE)
		tab_dist[tab_dist==0]=NA
		plot(seq_r, tab_dist, axes = FALSE, bty = "n", xlab = "", ylab = "",log="x",t="p",pch="*",ylim=c(0,2),cex=4)
		axis(side=4, at = pretty(c(0,2)),cex.axis=1.5)
		mtext("dist obs.", side=4, line=3,cex=1.)
	}

	if(sp==2){
		legend("topright",c(paste(expression(delta), "=10^",c("-3.2","-4","-5"),sep=""),"Sim.","Theory"),col=c("orchid","purple","orange","black","grey"),pch=c(16,17,18,NA,NA),lty=c(NA,NA,NA,1,2),bty="n",cex=1.5)
	}

}
dev.off()
