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

#Theoretical value of G
G_theoretical=function(gamma,r,lambda,Delta,C_0,Tmax=NA){ #U=0 in the absence of advection
  tau=0.0002
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

pcf_tot=c()
r_tot=c()
dominance_tot=c()
sp1_tot=c()
sp2_tot=c()

pdf("pcf_range_BBM_10m5.pdf")
for (nb_simu in 20:23){
f=read.table(paste("../simulation/pcf_",nb_simu,".txt",sep=""),sep=";",header=F,dec=".")
colnames(f)=c("r","sp1","sp2","pcf","dominance")
unique_sp=unique(f$sp1)
if (length(unique_sp)>length(colo)){
        colo=rainbow(length(unique_sp))
}
pcf_tot=c(pcf_tot,f$pcf[f$sp1==0&f$sp2==0])
r_tot=c(r_tot,f$r[f$sp1==0&f$sp2==0])
dominance_tot=c(dominance_tot,f$dominance[f$sp1==0&f$sp2==0])

#Count species
f_count=read.table(paste("../simulation/nb_indiv_",nb_simu,".txt",sep=""),header=F,sep=";",dec='.')
colnames(f_count)=c("species","abundance")

#Param simu
f_param=read.table(paste("../simulation/param_",nb_simu,".txt",sep=""),header=F,sep="=",dec='.')
f_param[,2]=as.numeric(as.character(f_param[,2]))

#####If we want to compare with spatstat in R
#  tab_repart=read.table(paste("../simulation/Spatial_distribution_7.txt",sep=""), header=F, sep=";")
#  colnames(tab_repart)=c("t","x","y","z","yfirst","ancestor","species")
#  tab_repart=subset(tab_repart,species==0)
#    ppp_repart=pp3(tab_repart$x,y=tab_repart$y,z=tab_repart$z,box3(c(0,1)))
#  pcf_poisson=pcf3est(ppp_repart,rmax=0.1,nrval=1000)
#       lines(pcf_poisson$r,pcf_poisson$trans,col="orchid")


delta_val=format(f_param[f_param[,1]=="delta",2],digits=2)
lambda=f_param[f_param[,1]=="growth_rate",2]
Delta=f_param[f_param[,1]=="Delta",2]
Utot=f_param[f_param[,1]=="Utot",2]
volume=f_param[f_param[,1]=="volume",2]

}

if(Utot==0.5){
        gamma=0.5
}else{
        stop("input gamma value")
}

	plot(1,1,log="xy",t="n",xlim=range(r_tot),ylim=range(pcf_tot[pcf_tot>0]),main=bquote(delta ~ "=" ~ .(delta_val)),axes=F)
        s1=1
        th_bbm=G_theoretical(gamma,unique(r_tot),lambda,Delta,f_count$abundance[f_count$species==unique_sp[s1]]/volume)
        if(nb_simu==2 || nb_simu==4){
                yl="g(r)"
        }else{
                yl=""
        }
        if(nb_simu<4){
                xl=""
        }else{
                xl="r"
        }
        if(nb_simu==2){
                m=bquote(delta ~ "=" ~ .(delta_val)~ "spatstat")
        }else{
                m=bquote(delta ~ "=" ~ .(delta_val))
        }
        lines(unique(r_tot),th_bbm,col="black",lty=2,lwd=2)
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()

        for(s2 in 1:length(unique_sp)){
                points(r_tot,pcf_tot,lty=1,col=colo[s2])
        }

dev.off()
