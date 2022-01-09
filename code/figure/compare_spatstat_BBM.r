rm(list=ls())
graphics.off()

library(spatstat)

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

BBM_cdf=function(r,gamma,Delta,Conc,lambda,t=NA){
  	tau=0.0002
	D=(Delta^2)/(2*tau)
        res=NA
        A=4/3*pi*r^3
        if(gamma>0){
                B=r^2/(6*D)
                C=(sqrt(gamma)/(6*sqrt(2)*(D^1.5)))*(r^3)*atan(r*sqrt(gamma/(2*D)))
                Dd=log((gamma*r^2)/(2*D)+1)/(6*gamma)
                E=(sqrt(gamma)*pi*r^3)/(12*sqrt(2)*D^1.5)
                res=A+(2*lambda/Conc)*(B+C+Dd-E)
        }else{
                B=r^2/2
                C=1/2*erf(r/(sqrt(8*D*t)))*(r^2-4*D*t)
                Dd=sqrt(2*D*t)*r/sqrt(pi)*exp(-r^2/(8*D*t))
                res=A+(lambda/(Conc*D))*(B-C-Dd)
        }
        return(res)
}

G_theoretical=function(gamma,r,lambda,Delta,C_0,Tmax=NA){ #U=0 in the absence of advection
  tau=0.0002
  D=(Delta^2)/(2*tau)
  if(gamma>0){
    tmp=lambda/(2*pi*C_0)*(sqrt(gamma)*atan(sqrt(gamma)*r/sqrt(2*D))/(2^1.5*D^1.5)+1/(2*D*r)-pi*sqrt(gamma)/(2^2.5*D^1.5))+1
  }else{ #This corresponds to the case U=0
        if(is.na(Tmax)){
                print("Tmax should have a value")
                stop()
        }
    tmp=lambda/(4*D*pi*C_0*r)*(1-erf(r/(sqrt(8*Tmax*D))))+1
  }
  return(tmp)
}



#####################################

colo=c("red","blue","grey","orange","violet")

nb_simu=80

#Parameters Diatoms
a_gamma=0.5

pdf("compare_BBM_spatstat_large_bandwidth.pdf",width=7.5)
par(mfrow=c(3,2),mar=c(4,4.5,1,4))

f=read.table(paste("../simulation/lambda_K_",nb_simu,"bis.txt",sep=""),sep=";",header=F,dec=".")
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

#Spatial distribution
tab_repart=read.table(paste("../simulation/Spatial_distribution_",nb_simu,".txt",sep=""), header=F, sep=";")
colnames(tab_repart)=c("t","x","y","z","yfirst","ancestor","species")

Concentration=f_count$abundance/f_param[f_param[,1]=="volume",2]

for (sp in unique(f_count$species)){
	if(sp==max(f_count$species)){
		xl="r"
	}else{
		xl=""
	}
  	tab_repart_tmp=subset(tab_repart,species==sp)
    	ppp_repart=pp3(tab_repart_tmp$x,y=tab_repart_tmp$y,z=tab_repart_tmp$z,box3(c(0,1)))

 	if(1==0){	
        dist=unique(c(pairdist(ppp_repart)))

        a=sort(dist[dist>0],decreasing=F)
        a=a[1:1000]
        f_plot=subset(f,sp1==sp&sp2==sp)
        seq_r=unique(f_plot$r)
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

	f_plot=subset(f,sp1==sp&sp2==sp)

#Ripley's function
  	K_spat=K3est(ppp_repart,rmax=max(f$r),nrval=1000)
	
	#Theory
	th_bbm=BBM_cdf(f_plot$r,a_gamma,a_Delta,Concentration[sp+1],a_lambda)

	plot(0.1,0.1,log="xy",xlab=xl,ylab="K",axes=F,xlim=range(f$r[f$r>0]),ylim=range(c(range(th_bbm),range(f$K[f$K>0]))),cex.lab=1.5,t="n")
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()
	points(K_spat$r,K_spat$trans,col="orchid",cex=1.5,pch=16)
	lines(f_plot$r,f_plot$K,col=colo[sp+1])
	lines(f_plot$r,th_bbm,lty=2,col="grey")
	print(mean(th_bbm))
	if(sp==0){
		legend("bottomright",c("Spatstat","Simulations","Theory"),lty=c(NA,1,2),pch=c(16,NA,NA),col=c("orchid",colo[sp+1],"grey"),bty="n")
	}
#PCF
  	pcf_spat=pcf3est(ppp_repart,rmax=max(f$r),nrval=1000,delta=10^(-5))

	#Theory
	th_bbm=G_theoretical(a_gamma,f_plot$r,a_lambda,a_Delta,Concentration[sp+1])

	plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab="PCF",axes=F,xlim=range(f$r[f$r>0]),cex.lab=1.5,ylim=range(f$pcf[f$pcf>0]))
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()
	points(pcf_spat$r,pcf_spat$trans,col="orchid",cex=1.5,pch=16)
        lines(f_plot$r,f_plot$pcf,col=colo[sp+1])
	lines(f_plot$r,th_bbm,lty=2,col="grey")

	if(1==0){
	par(new = TRUE)
        tab_dist[tab_dist==0]=NA
        plot(seq_r, tab_dist, axes = FALSE, bty = "n", xlab = "", ylab = "",log="x",t="p",pch="*",ylim=c(0,2),cex=3)
        axis(side=4, at = pretty(c(0,2)),cex.axis=1.5)
        mtext("dist obs.", side=4, line=3,cex=1.)
	}

}
dev.off()
