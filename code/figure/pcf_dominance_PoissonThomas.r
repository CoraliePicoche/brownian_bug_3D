rm(list=ls())
graphics.off()

############################### Utilitary functions

#Ticks for log scale from https://stackoverflow.com/questions/47890742/logarithmic-scale-plot-in-r
log10Tck <- function(side, type){
  lim <- switch(side,
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceiling(lim[2])
  return(switch(type,
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}

#####################################

colo=c("red","blue","grey","orange","violet")

pdf("metrics_PoissonThomas.pdf",width=18.5,height=9)
par(mfrow=c(2,4),cex=1.2,mar=c(5,5,0.5,1.5))

#Poisson Simu
nb_simu_tot=1:4

f_tot=matrix(,0,5)
colnames(f_tot)=c("r","sp1","sp2","pcf","dominance")
f_count_tot=matrix(,0,2)
colnames(f_count_tot)=c("species","abundance")
f_param_tot=matrix(,0,2)
colnames(f_param_tot)=c("name","value")

for (nb_simu in nb_simu_tot){
	f=read.table(paste("../simulation/pcf_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
	colnames(f)=c("r","sp1","sp2","pcf","dominance")
	unique_sp=unique(f$sp1)
	if (length(unique_sp)>length(colo)){
		colo=rainbow(length(unique_sp))
	}
	f_tot=rbind(f_tot,f)

	#Count species
	f_count=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";",dec='.')
	colnames(f_count)=c("species","abundance")
	f_count_tot=rbind(f_count_tot,f_count)

	#Param simu
	f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=",dec='.')
	f_param[,2]=as.numeric(as.character(f_param[,2]))
	colnames(f_param)=c("name","value")
	f_param_tot=rbind(f_param_tot,f_param)

} #end loop on nb_simu

for (s1 in 1:length(unique_sp)){
	if(s1==1){
		yl="g(r)"
	}else{
		yl=""
	}	
	plot(0.1,0.1,t="n",log="xy",xlab="",ylab=yl,axes=F,xlim=range(f_tot$r[f_tot$r>0]),cex.lab=1.5,ylim=range(f_tot$pcf[f_tot$pcf>0],na.rm=T))
	if(s1==2){
		legend("topleft",c(paste("SP=",unique_sp),"theory"),col=c(colo[1:length(unique_sp)],"black"),lty=c(rep(1,length(unique_sp)),2),lwd=2,bty="n",pch=c(rep(1,length(unique_sp),NA)))
	}
	f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])
  	lines(f_plot$r,rep(1,length(f_plot$r)),col="black",lty=2,lwd=2)
	axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
	box()

	for(s2 in 1:length(unique_sp)){
		f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s2])
		lines(f_plot$r,f_plot$pcf,lty=1,col=colo[s1],lwd=2)
		points(f_plot$r,f_plot$pcf,lty=1,col=colo[s2])
	}
}

plot(1,1,t="n",xlim=range(f_tot$r[f$r>0]),ylim=range(f_tot$dominance,na.rm=T),xlab="",ylab="dominance",log="x",cex.lab=1.5,cex.axis=1.5,axes=F)
axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2,tcl=0.5,cex.axis=1.5) # left
box()

for (s1 in 1:length(unique_sp)){
	f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])
	points(f_plot$r,f_plot$dominance,col=colo[s1])
}

#Thomas distribution
nb_simu_tot=5:8 #These are the ones with the delta spatstat

f_tot=matrix(,0,5)
colnames(f_tot)=c("r","sp1","sp2","pcf","dominance")
f_count_tot=matrix(,0,2)
colnames(f_count_tot)=c("species","abundance")
f_param_tot=matrix(,0,2)
colnames(f_param_tot)=c("name","value")

for (nb_simu in nb_simu_tot){
        f=read.table(paste("../simulation/pcf_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
        colnames(f)=c("r","sp1","sp2","pcf","dominance")
        unique_sp=unique(f$sp1)
        if (length(unique_sp)>length(colo)){
                colo=rainbow(length(unique_sp))
        }
        f_tot=rbind(f_tot,f)

        #Count species
        f_count=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";",dec='.')
        colnames(f_count)=c("species","abundance")
        f_count_tot=rbind(f_count_tot,f_count)

        #Param simu
        f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=",dec='.')
        f_param[,2]=as.numeric(as.character(f_param[,2]))
        colnames(f_param)=c("name","value")
        f_param_tot=rbind(f_param_tot,f_param)

} #end loop on nb_simu


sigma=f_param_tot[f_param_tot$name=="sigma","value"]
if(length(unique(sigma))>1){
	stop("sigma were different")
}
sigma=unique(sigma)

N_parent=f_param_tot[f_param_tot$name=="N_parent","value"]
if(length(unique(N_parent))>1){
	stop("N_parent were different")
}
N_parent=unique(N_parent)

volume=f_param_tot[f_param_tot$name=="volume","value"]
if(length(unique(volume))>1){
	stop("volume were different")
}
volume=unique(volume)
C_parent=N_parent/volume
th_thomas=1+exp(-unique(f_tot$r)^2/(4*sigma^2))*(4*pi*sigma^2)^(-3/2)*1/C_parent

for (s1 in 1:length(unique_sp)){
        if(s1==1){
                yl="g(r)"
        }else{
                yl=""
        }
        plot(0.1,0.1,t="n",log="xy",xlab="r",ylab=yl,axes=F,xlim=range(f_tot$r[f_tot$r>0]),cex.lab=1.5,ylim=range(c(f_tot$pcf[f_tot$pcf>0],th_thomas),na.rm=T))
        f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()

        for(s2 in 1:length(unique_sp)){
                f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s2])
                lines(f_plot$r,f_plot$pcf,lty=1,col=colo[s1],lwd=2)
                points(f_plot$r,f_plot$pcf,lty=1,col=colo[s2])
        }
	lines(unique(f_tot$r),th_thomas,col="black",lty=2,lwd=2)
}

plot(1,1,t="n",xlim=range(f_tot$r[f$r>0]),ylim=range(f_tot$dominance,na.rm=T),xlab="r",ylab="dominance",log="x",cex.lab=1.5,cex.axis=1.5,axes=F)
axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(2,tcl=0.5,cex.axis=1.5) # left
box()
for (s1 in 1:length(unique_sp)){
        f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])
        points(f_plot$r,f_plot$dominance,col=colo[s1])
}

dev.off()
