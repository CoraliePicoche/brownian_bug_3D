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
a_pch=c(16,17)
a_cex=c(1.25,0.75)

nb_simu=80

a_gamma=0.5

pdf("BBM_interspecies.pdf",width=10,height=7)
par(mfcol=c(2,3),mar=c(4,5,3,4))

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

unique_sp=unique(f_count$species)

#tab_repart_tmp=subset(tab_repart,species==sp)
ppp_repart=pp3(tab_repart$x,y=tab_repart$y,z=tab_repart$z,box3(c(0,1)))

for (sp in unique_sp){
	osp=setdiff(unique_sp,sp)
	xl="r"
	
	#Ripley's function
	if(sp==0){
		yl="K"
	}else{
		yl=""
	}

	plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab=yl,axes=F,xlim=range(f$r[f$r>0]),cex.lab=1.5,ylim=range(f$K[f$K>0]),main=paste("Species",sp))
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()

	p=0
	colo_leg=c()	
	for(sp2bis in osp){
	p=p+1
	colo_leg=c(colo_leg,colo[sp2bis+1])
	f_plot=subset(f,sp1==sp&sp2==sp2bis)
	seq_r=unique(f_plot$r)
	
	#Theory
	#th_bbm_mono=BBM_cdf(f_plot$r,a_gamma,a_Delta,Concentration[sp+1],a_lambda)
	th_bbm_inter=4/3*pi*(f_plot$r^3)

        points(f_plot$r,f_plot$K,col=colo[sp2bis+1],pch=a_pch[p],cex=a_cex[p])
	#lines(f_plot$r,th_bbm_mono,lty=2,col="grey")
	lines(f_plot$r,th_bbm_inter,lty=2,col="black",lwd=2)

	legend("topleft",c(paste("Sp",osp),"Theory"),col=c(colo_leg,"black"),pch=c(a_pch,NA),lty=c(NA,NA,2),bty="n",cex=1.5)

	}#end other species

	#PCF
        if(sp==0){
                yl="PCF"
        }else{
                yl=""
        }

        plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab=yl,axes=F,xlim=range(f$r[f$r>0]),cex.lab=1.5,ylim=range(f$pcf[f$pcf>0]))
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()

        p=0
        colo_leg=c()
        for(sp2bis in osp){
        p=p+1
        colo_leg=c(colo_leg,colo[sp2bis+1])
        f_plot=subset(f,sp1==sp&sp2==sp2bis)
        seq_r=unique(f_plot$r)

        #Theory
	#th_bbm=G_theoretical(a_gamma,f_plot$r,a_lambda,a_Delta,Concentration[sp+1])
        th_bbm_inter=rep(1,length(f_plot$r))

        points(f_plot$r,f_plot$pcf,col=colo[sp2bis+1],pch=a_pch[p])
        #lines(f_plot$r,th_bbm_mono,lty=2,col="grey")
        lines(f_plot$r,th_bbm_inter,lty=2,col="black")

        }#end other species


	} #end sp
dev.off()
