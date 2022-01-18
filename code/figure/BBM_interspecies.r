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

colo=c("red","blue","black")
a_pch=c(16,17)
a_cex=c(1.25,0.75)

nb_simu_list=c(150,151) 
a_main_list=c("Advection","No advection")

pdf("BBM_interspecies_micro.pdf",width=10,height=7)
par(mfcol=c(3,2),mar=c(4,5,3,4))

for(i in 1:length(nb_simu_list)){
	    nb_simu=nb_simu_list[i]

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

#Spatial distribution
tab_repart=read.table(paste("../simulation/Spatial_distribution_",nb_simu,".txt",sep=""), header=F, sep=";")
colnames(tab_repart)=c("t","x","y","z","yfirst","ancestor","species")

Concentration=f_count$abundance/f_param[f_param[,1]=="volume",2]

unique_sp=unique(f_count$species)

for (sp in unique_sp){
	osp=setdiff(unique_sp,sp)
	xl="r"
	
	#Ripley's function
	if(sp==0){
		yl="K"
		ml=paste(a_main_list[i],"\n Sp =",sp)
	}else{
		yl=""
		ml=paste("Sp =",sp)
	}


	#Theory
	seq_r=unique(f$r)
	th_bbm_inter=4/3*pi*(seq_r^3)
	
	plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab=yl,axes=F,xlim=range(f$r[f$r>0]),cex.lab=1.5,ylim=range(c(f$K[f$K>0]),th_bbm_inter),main=ml)
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

        points(f_plot$r,f_plot$K,col=colo[sp2bis+1],pch=a_pch[p],cex=a_cex[p])
	lines(seq_r,th_bbm_inter,lty=2,col="grey",lwd=2)

	if(sp==0){
		legend("topleft",c(paste("Sp",osp),"Theory"),col=c(colo_leg,"grey"),pch=c(a_pch,NA),lty=c(NA,NA,2),bty="n",cex=1.5)
		mtext("Micro", outer = TRUE, cex = 1.5,line=-2.)
	}

	}#end other species

	} #end sp
		}#end nb_simu
dev.off()
