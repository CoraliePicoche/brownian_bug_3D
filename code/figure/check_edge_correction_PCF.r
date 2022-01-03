rm(list=ls())
graphics.off()

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
thomas_cdf=function(r,sigma){
        a=2/(sigma*sqrt(pi))*(erf(r/(2*sigma))-r/(sigma*sqrt(pi))*exp(-(r^2/(4*sigma^2))))
        return(a)
}

#####################################

colo=c("red","blue","grey")

pdf("PCF_Poisson_edge_correction.pdf",width=10,height=10)
par(mfrow=c(3,3),mar=c(4,4.5,3,1))

#sim=list(101,103,105)#Thomas 
sim=list(100,102,104)#Poisson
type=c("Translation","Torus","Both") 

for (species in 1:3){
	s1=species
	ss=s1
	if(s1==3){
		xl="r"
	}else{
		xl=""
	}
for (s in 1:length(sim)){
	print(type[s])
	nb_simu_tot=sim[s][[1]]
	if(s==1){
		yl="PCF"
	}else{
		yl=""
	}

	f_tot=matrix(,0,7)
	colnames(f_tot)=c("r","sp1","sp2","pcf","dominance","lambda_K","K")
	f_count_tot=matrix(,0,2)
	colnames(f_count_tot)=c("species","abundance")
	f_param_tot=matrix(,0,2)
	colnames(f_param_tot)=c("name","value")

	for (nb_simu in nb_simu_tot){
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
	th_thomas=1+exp(-unique(f_tot$r)^2/(4*sigma^2))*(4*pi*sigma^2)^(-3/2)*1/(N_parent/volume)
	th_poisson=rep(1,length(unique(f_tot$r)))

	plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab=yl,axes=F,xlim=c(10^(-4),0.1),cex.lab=1.5,ylim=range(f_tot$pcf)+10^(-4),main=paste("Species=",s1," correct=",type[s],sep=""))
	#plot(0.1,0.1,t="n",log="xy",xlab=xl,ylab=yl,axes=F,xlim=c(10^(-4),0.1),cex.lab=1.5,ylim=range(th_poisson),main=paste("Species=",s1," correct=",type[s],sep=""))
	axis(1, at=log10Tck('x','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
	axis(2, at=log10Tck('y','major'), tcl= 0.5,cex.axis=1.5) # bottom
	axis(2, at=log10Tck('y','minor'), tcl= 0.1, labels=NA) # bottom
	box()

	f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s1])

	points(f_plot$r,f_plot$pcf,pch=16,col=colo[s1])
	colors=colo[s1]

	print(mean(f_plot$pcf))
	for(s2 in 1:length(unique_sp)){
		if(s1!=s2){
			colors=c(colors,colo[s2])
			if(s==2){ #Only for the legend
				ss=c(ss,s2)
			}
			f_plot=subset(f_tot,sp1==unique_sp[s1]&sp2==unique_sp[s2])
			lines(f_plot$r,f_plot$pcf,lty=1,col=colo[s2])
		}
	}
	if(s==2){	
		#legend("topleft",c(paste("S=",s1," x S=",ss,sep=""),"Theory mono-specific", "Theory cross-specific"),col=c(colors,"black","black"),bty="n",pch=c(16,NA,NA,NA,NA),lty=c(NA,1,1,2,3),lwd=2)
		legend("topleft",c(paste("S=",s1," x S=",ss,sep=""),"Theory"),col=c(colors,"black"),bty="n",pch=c(16,NA,NA,NA),lty=c(NA,1,1,2),lwd=2)
	}
	lines(unique(f_tot$r),th_poisson,lty=3,lwd=2)
	#lines(unique(f_tot$r),th_thomas+th_poisson,lty=2,lwd=2)
}
}
dev.off()
