rm(list=ls())
graphics.off()

library(spatstat)
library(pracma)

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

BBM_cdf=function(r,gamma,Delta,Conc,lambda,tau,t=NA){
	D=(Delta^2)/(2*tau)
        res=NA
        A=4/3*pi*r^3
	t=t*tau
        if(gamma>0){
                B=r^2/(6*D)
                C=(sqrt(gamma)/(6*sqrt(2)*(D^1.5)))*(r^3)*atan(r*sqrt(gamma/(2*D)))
                Dd=log((gamma*r^2)/(2*D)+1)/(6*gamma)
                E=(sqrt(gamma)*pi*r^3)/(12*sqrt(2)*D^1.5)
                res=A+(2*lambda/Conc)*(B+C+Dd-E)
        }else{
                B=r^2/2
                C=1/2*erf(r/(sqrt(8*D*t)))*(r^2-4*D*t)
		print(D)
		print(head(erf(r/(sqrt(8*D*t)))))
                Dd=sqrt(2*D*t)*r/sqrt(pi)*exp(-r^2/(8*D*t))
                res=A+(lambda/(Conc*D))*(B-C-Dd)
        }
        return(res)
}

#####################################

colo=c("red","blue","black")

nb_simu_list=c(150,151)

a_gamma_list=c(1231,0)

a_main_list=c("Advection","No advection")

pdf("K_micro.pdf",width=7.5)
par(mfcol=c(3,2),mar=c(4,4.5,5,4))


for(i in 1:length(nb_simu_list)){

	    nb_simu=nb_simu_list[i]
	    a_gamma=a_gamma_list[i]

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
a_lambda=f_param[f_param[,1]=="growth_rate",2]
a_Delta=f_param[f_param[,1]=="Delta",2]
a_tau=f_param[f_param[,1]=="tau",2]

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

        dist=unique(c(pairdist(ppp_repart)))

	if(1==0){
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
	#Theory
	th_bbm=BBM_cdf(f_plot$r,a_gamma,a_Delta,Concentration[sp+1],a_lambda,a_tau,t=1000)
  	
#	K_spat=K3est(ppp_repart,rmax=max(f$r),nrval=500)

	if(sp==0){
		ml=paste(a_main_list[i],"\n Sp =",sp)
	}else{
		ml=paste("Sp =",sp)
	}

	plot(0.1,0.1,log="xy",xlab=xl,ylab="K",axes=F,xlim=range(f$r[f$r>0]),ylim=range(c(range(th_bbm),range(f$K[f$K>0]))),cex.lab=1.5,t="n",main=ml)
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()

	if(sp==0){
		mtext("Micro", outer = TRUE, cex = 1.5,line=-2.)
		legend("bottomright",c("Simulations","Theory","10 radii"),lty=c(1,2,3),pch=c(NA,NA),col=c(colo[sp+1],"grey","grey"),bty="n",lwd=2)
	}

	#	points(K_spat$r,K_spat$trans,col="orchid",cex=1.5,pch=16)
	lines(f_plot$r,f_plot$K,col=colo[sp+1],lwd=2)
	lines(f_plot$r,th_bbm,lty=2,col="grey",lwd=2)
	print(mean(th_bbm))

	if(1==0){
	par(new = TRUE)
        tab_dist[tab_dist==0]=NA
        plot(seq_r, tab_dist, axes = FALSE, bty = "n", xlab = "", ylab = "",log="x",t="p",pch="*",ylim=c(0,2),cex=3)
        axis(side=4, at = pretty(c(0,2)),cex.axis=1.5)
        mtext("dist obs.", side=4, line=3,cex=1.)
	}

#	abline(v=1.5*10^(-4)*10,col="grey",lwd=2,lty=3)
	abline(v=25*10^(-4)*10,col="grey",lwd=2,lty=3)

}
}
dev.off()
