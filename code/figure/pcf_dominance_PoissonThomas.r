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

#Theoretical value of G
G_theoretical=function(gamma,r,lambda,Delta,C_0,Tmax=NA){ #U=0 in the absence of advection
  tau=1
  D=Delta^2/(2*tau)
  if(gamma>0){
    tmp=lambda/(4*pi*C_0)*(sqrt(gamma)*atan(sqrt(gamma)*r/sqrt(2*D))/(2^1.5*D^1.5)+1/(2*D*r)-pi*sqrt(gamma)/(2^2.5*D^1.5))+1
  }else{ #This corresponds to the case U=0
        if(is.na(Tmax)){
                print("Tmax should have a value")
                stop()
        }
    tmp=2*lambda/(8*D*pi*C_0*r)*(1-erf(r/(2^1.5*sqrt(Tmax*D))))+1
  }
  return(tmp)
}

#####################################

colo=c("red","blue","grey","orange","violet")

pdf("metrics_PoissonThomas.pdf",width=15,height=9)
par(mfrow=c(2,4),cex=1.,mar=c(5,5,0.5,1))

nb_simu=1 #Poisson
f=read.table(paste("../simulation/pcf_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f)=c("r","sp1","sp2","pcf","dominance")
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

dd=c()
for (s1 in 1:length(unique_sp)){
	if(s1==1){
		yl="g(r)"
	}else{
		yl=""
	}	
	plot(0.1,0.1,t="n",log="xy",xlab="",ylab=yl,axes=F,xlim=range(f$r[f$r>0]),cex.lab=1.5,ylim=range(f$pcf))
	f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s1])
  	lines(f_plot$r,rep(1,length(f_plot$r)),col="black",lty=2,lwd=2)
	axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	axis(2,tcl=0.5,cex.axis=1.5) # left
	box()

	for(s2 in 1:length(unique_sp)){
		f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s2])
		lines(f_plot$r,f_plot$pcf,lty=1,col=colo[s1],lwd=2)
		points(f_plot$r,f_plot$pcf,lty=1,col=colo[s2])
	}
}

plot(1,1,t="n",xlim=range(f$r[f$r>0]),ylim=range(f$dominance,na.rm=T),xlab="r",ylab="dominance",log="x",cex.lab=1.5,cex.axis=1.5)
for (s1 in 1:length(unique_sp)){
	f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s1])
	points(f_plot$r,f_plot$dominance,col=colo[s1])
}
legend("topright",paste("SP=",unique_sp),col=colo,lty=1,bty="n")

nb_simu=2 #Thomas
f=read.table(paste("../simulation/pcf_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f)=c("r","sp1","sp2","pcf","dominance")
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

C_parent=f_param[f_param[,1]=="N_parent",2]/f_param[f_param[,1]=="volume",2]
sigma=f_param[f_param[,1]=="sigma",2]

dd=c()
for (s1 in 1:length(unique_sp)){
	if(s1==1){
		yl="g(r)"
	}else{
		yl=""
	}	
	plot(0.1,0.1,t="n",log="xy",xlab="r",ylab=yl,axes=F,xlim=range(f$r[f$r>0]),cex.lab=1.5,ylim=range(f$pcf))
        f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s1])
  	th_thomas=1+exp(-f_plot$r^2/(4*sigma^2))*(4*pi*sigma^2)^(-3/2)*1/C_parent
        lines(f_plot$r,th_thomas,col="black",lty=2,lwd=2)
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2,tcl=0.5,cex.axis=1.5) # left
        box()

        for(s2 in 1:length(unique_sp)){
                f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s2])
                lines(f_plot$r,f_plot$pcf,lty=1,col=colo[s1],lwd=2)
                points(f_plot$r,f_plot$pcf,lty=1,col=colo[s2])
        }
}

plot(1,1,t="n",xlim=range(f$r[f$r>0]),ylim=range(f$dominance,na.rm=T),xlab="r",ylab="dominance",log="x",cex.axis=1.5,cex.lab=1.5)
for (s1 in 1:length(unique_sp)){
        f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s1])
        points(f_plot$r,f_plot$dominance,col=colo[s1])
}
legend("topright",paste("SP=",unique_sp),col=colo,lty=1,bty="n")

dev.off()
