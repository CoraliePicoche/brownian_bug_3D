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

colo=c("red","blue","green","orange","grey")
volume=10

pdf("cross_pcf_2sp_BBM_test_pcf.pdf",width=16,height=12)
par(mfrow=c(2,3))

f=read.table("../simulation/pcf_BBM_kernel_2sp_realistic_values_20000_U0p5_test_pcf.txt",sep=";",header=F,dec=".")
colnames(f)=c("r","sp1","sp2","pcf","dominance")
unique_sp=unique(f$sp1)
if (length(unique_sp)>length(colo)){
	colo=rainbow(length(unique_sp))
}

#Count species
f_count=read.table("../simulation/nb_indiv_BBM_kernel_2sp_realistic_values_20000_U0p5_test_pcf.txt",header=F,sep=";",dec='.')
colnames(f_count)=c("species","abundance")

dd=c()

for (s1 in 1:length(unique_sp)){
	plot(0.1,0.1,t="n",log="xy",xlab="r",ylab="",axes=F,ylim=c(0.5,max(f$pcf,na.rm=T)),xlim=range(f$r)+10^(-4),cex.lab=1.5)
	f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s1])
  	lines(f_plot$r,G_theoretical(0.5,f_plot$r,0.0002,7*10^(-5),f_count$abundance[f_count$species==(s1-1)]/volume),col="black",lty=2,lwd=2)
	axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	axis(2,tcl=0.5,cex.axis=1.5) # left
	box()

	for(s2 in 1:length(unique_sp)){
		f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s2])
		if(s1==s2){
		lines(f_plot$r,f_plot$pcf,lty=1,col=colo[s1],lwd=2)
		}else{
			dd=c(dd,f_plot$pcf[f_plot$r<0.01])
		}
		points(f_plot$r,f_plot$pcf,lty=1,col=colo[s2])
	}
}

plot(1,1,t="n",xlim=range(f$r)+0.00001,ylim=range(f$dominance,na.rm=T),xlab="r",ylab="dominance",log="x")
for (s1 in 1:length(unique_sp)){
	f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s1])
	points(f_plot$r,f_plot$dominance,col=colo[s1])
}
legend("topright",paste("SP=",unique_sp),col=colo,lty=1,bty="n")

#hist(log(dd),main="log(pcf),r<0.01",xlab="log(pcf)")
#abline(v=0)

dev.off()
