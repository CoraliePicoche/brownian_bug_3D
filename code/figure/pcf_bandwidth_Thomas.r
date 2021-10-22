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

#####################################

colo=c("red","blue","grey","orange","violet")

pdf("pcf_bandwidth_Thomas.pdf")
par(mfrow=c(2,3),mar=c(4,4,3,0.5))

for (nb_simu in c(2,7,3,4,5,6)){
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

#####If we want to compare with spatstat in R
#  tab_repart=read.table(paste("../simulation/Spatial_distribution_7.txt",sep=""), header=F, sep=";")
#  colnames(tab_repart)=c("t","x","y","z","yfirst","ancestor","species")
#  tab_repart=subset(tab_repart,species==0)
#    ppp_repart=pp3(tab_repart$x,y=tab_repart$y,z=tab_repart$z,box3(c(0,1)))
#  pcf_poisson=pcf3est(ppp_repart,rmax=0.1,nrval=1000)
#	lines(pcf_poisson$r,pcf_poisson$trans,col="orchid")


C_parent=f_param[f_param[,1]=="N_parent",2]/f_param[f_param[,1]=="volume",2]
sigma=f_param[f_param[,1]=="sigma",2]
delta_val=format(f_param[f_param[,1]=="delta",2],digits=2)

  	th_thomas=1+exp(-unique(f$r)^2/(4*sigma^2))*(4*pi*sigma^2)^(-3/2)*1/C_parent
s1=1
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
	plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab=yl,axes=F,xlim=range(f$r[f$r>0]),cex.lab=1.5,ylim=c(1,max(c(th_thomas),na.rm=T)),main=m)
        f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s1])
        lines(f_plot$r,th_thomas,col="black",lty=2,lwd=2)
        axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
        axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
        axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
        box()

        for(s2 in 1:length(unique_sp)){
                f_plot=subset(f,sp1==unique_sp[s1]&sp2==unique_sp[s2])
                lines(f_plot$r,f_plot$pcf,lty=1,col=colo[s1],lwd=2)
                points(f_plot$r,f_plot$pcf,lty=1,col=colo[s2])
        }

}
dev.off()
