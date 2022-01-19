rm(list=ls())
graphics.off()

library("spatstat")
source("utilitary_functions.r")
source("theoretical_functions.r")


colo=c("red","blue","black")

pdf("K_PCF_Thomas_edge_correction_large_bandwidth.pdf",width=8,height=12)
par(mfrow=c(3,2))

nb_simu=107

for (species in 1:3){
	s1=species
	ss=s1
	if(s1==3){
		xl="r"
	}else{
		xl=""
	}
	

	f_tot=matrix(,0,7)
	colnames(f_tot)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
	f_count_tot=matrix(,0,2)
	colnames(f_count_tot)=c("species","abundance")
	f_param_tot=matrix(,0,2)
	colnames(f_param_tot)=c("name","value")

		f=read.table(paste("../simulation/lambda_K_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
		colnames(f)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
		f_tot=rbind(f_tot,f)
		unique_sp=unique(f_tot$sp1)

		#Count species
		f_count=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";",dec='.')
		colnames(f_count)=c("species","abundance")
		f_count_tot=rbind(f_count_tot,f_count)

		#Param simu
		f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=",dec='.')
		f_param[,2]=as.numeric(as.character(f_param[,2]))
		colnames(f_param)=c("name","value")
		f_param_tot=rbind(f_param_tot,f_param)
	sigma=f_param_tot[f_param_tot$name=="sigma","value"]
	if(length(unique(sigma))>1){
		stop("sigma were different")
	}
	sigma=unique(sigma)

	N_parent=f_param_tot[f_param_tot$name=="N_parent","value"]
	print(N_parent)
	if(length(unique(N_parent))>1){
		stop("N_parent were different")
	}
	N_parent=unique(N_parent)

	volume=f_param_tot[f_param_tot$name=="volume","value"]
	if(length(unique(volume))>1){
		stop("volume were different")
	}
	volume=unique(volume)
	th_thomas=thomas_cdf(unique(f_tot$r),sigma)/(N_parent/volume)
	th_poisson=4/3*pi*unique(f_tot$r)^3

	plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab="K",axes=F,xlim=c(10^(-4),0.1),cex.lab=1.5,ylim=range(f_tot$K[f_tot$K>0]),main=paste("Species=",unique_sp[s1],sep=""))
	axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
	box()

	f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])

#Spatstat
	tab_repart=read.table(paste("../simulation/Spatial_distribution_",nb_simu,".txt",sep=""), header=F, sep=";")
	colnames(tab_repart)=c("t","x","y","z","yfirst","ancestor","species")

	tab_repart_tmp=subset(tab_repart,species==unique_sp[s1])
	ppp_repart=pp3(tab_repart_tmp$x,y=tab_repart_tmp$y,z=tab_repart_tmp$z,box3(c(0,1)))
  	K_spat=K3est(ppp_repart,rmax=max(f$r),nrval=1000)
	points(K_spat$r,K_spat$trans,col="orchid",pch=16,cex=1.5)

	points(f_plot$r,f_plot$K,pch=16,col=colo[s1],cex=0.5)
	colors=colo[s1]

	for(s2 in 1:length(unique_sp)){
		if(s1!=s2){
			colors=c(colors,colo[s2])
				ss=c(ss,unique_sp[s2])
			f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s2])
			print(paste("SP",s1," X ",s2,sep=""))
			print(sum(which(f_plot$K>0)))
			lines(f_plot$r,f_plot$K,lty=1,col=colo[s2])
		}
	}
		legend("topleft",c(paste("S=",s1," x S=",ss,sep=""),"Theory: intraspecific", "Theory: interspecific","Spatstat"),col=c(colors,"grey","grey","orchid"),bty="n",pch=c(16,NA,NA,NA,NA,16),lty=c(NA,1,1,2,3,NA),lwd=c(rep(2,4),3,2),cex=1.25)
	lines(unique(f_tot$r),th_poisson,lty=3,lwd=3,col="grey")
	lines(unique(f_tot$r),th_thomas+th_poisson,lty=2,lwd=2,col="grey")

        th_thomas=1+exp(-unique(f_tot$r)^2/(4*sigma^2))*(4*pi*sigma^2)^(-3/2)*1/(N_parent/volume)
        th_poisson=rep(1,length(unique(f_tot$r)))
  	
	PCF_spat=pcf3est(ppp_repart,rmax=max(f$r),nrval=1000,delta=10^(-5))

        plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab="PCF",axes=F,xlim=c(10^(-4),0.1),cex.lab=1.5,ylim=range(f_tot$pcf[f_tot$pcf>0]),main=paste("Species=",unique_sp[s1],sep=""))
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()

        f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])

	points(PCF_spat$r,PCF_spat$trans,col="orchid",pch=16,cex=1.5)
        points(f_plot$r,f_plot$pcf,pch=16,col=colo[s1],cex=0.5)
        colors=colo[s1]

        for(s2 in 1:length(unique_sp)){
                if(s1!=s2){
                        f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s2])
                        print(paste("SP",s1," X ",s2,sep=""))
                        print(sum(which(f_plot$pcf>0)))
                        lines(f_plot$r,f_plot$pcf,lty=1,col=colo[s2])
                }
        }
        lines(unique(f_tot$r),th_poisson,lty=3,lwd=3,col="grey")
        lines(unique(f_tot$r),th_thomas+th_poisson,lty=2,lwd=2,col="grey")

	#a=loess(f_plot$pcf~f_plot$r,span=1)
	#lines(f_plot$r,predict(a,newdata=f_plot$r),col="orchid",lwd=2)
}
dev.off()

